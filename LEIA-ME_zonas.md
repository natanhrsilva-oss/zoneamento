# Zonas de amostragem direcionada · TruBio

Ferramenta interna: recebe o limite do talhão e propõe um **grid variável** (uma zona = uma amostra composta),
com base no histórico de vigor (NDVI Sentinel-2) e no relevo (Copernicus DEM). Só fontes públicas e gratuitas.

## Uso
    pip install -r requirements_zonas.txt
    streamlit run app_zonas.py                      # aplicativo
    python zonas_core.py talhao.kml --kml zonas.kml # linha de comando (--camadas, --painel, --csv, --demo)
    python zonas_core.py --selftest                 # 46 verificações, sem rede

Arquivos: `app_zonas.py` (interface), `zonas_core.py` (motor), `zonas_persistencia.py` (funções de apoio e extração
por célula para validação), `zonas_selftest.py` (testes).

## Entrada e saída
- Entrada: KML, KMZ, shapefile (.zip ou .shp+.dbf+.shx+.prj) ou GeoJSON, em qualquer CRS (convertido para WGS84).
  Se o arquivo for um grid de amostragem já existente, marque "unir em um talhão por arquivo".
- Saída para a plataforma: **KML (WGS84) com um Placemark/Polygon por zona**, nome `{talhao}_Z{n:02d}` (configurável).
  O modo "KML limpo" grava só nome e geometria. Também: KML de camadas, CSV das zonas e PNG do painel.

## Vistas do aplicativo
1. RGB do talhão (uma cena Sentinel-2 limpa e recente, sem mapa-base externo).
2. Clusterização: as zonas propostas (cores e áreas).
3. Zonas de estabilidade: sem desvios relevantes sobre a **mediana** do talhão.
4. Zonas de interesse: queda de vigor, persistente ou recente.
No mapa interativo as quatro são camadas liguáveis. O painel traz as quatro lado a lado.

## Definições (ajustáveis na barra lateral)
- Série: janelas de safra V (out–mar) e I (abr–set), 2017 em diante. NDVI no percentil 90 por pixel e janela.
  Desvio relativo = (NDVI − mediana do talhão na janela) / mediana. Janelas sem lavoura (mediana < 0,5) são descartadas.
- **Estável:** |desvio médio| ≤ limiar/2 e desvio acima do limiar (+ ou −) em ≤ 15% das janelas.
- **Queda persistente:** desvio < −limiar em ≥ 50% das janelas. **Queda recente:** < −limiar em ≥ 2/3 das últimas 3
  janelas e normal antes. Limiar padrão: 10%.
- Zonas: N = área do talhão ÷ densidade-alvo. Limites de área **flexíveis** (penalizados, não rígidos) e compacidade
  ajustável. Zonas abaixo de 2 ha são fundidas (inviáveis de amostrar).
- **Suavização dos contornos** (padrão 0,8; 0 desliga, 1 é o máximo): em duas etapas, sempre sobre a rede de zonas inteira.
  (1) a partição é suavizada numa grade 4× mais fina; (2) as fronteiras compartilhadas são suavizadas como curvas, em arcos
  entre os pontos onde três zonas se encontram (esses pontos se acomodam ao centro das curvas e levam junto as pontas dos
  arcos vizinhos). Cada fronteira é suavizada uma só vez, então as zonas seguem encaixadas, sem frestas nem sobreposição, e a
  borda real do talhão não é alterada. Vale também para as camadas de estabilidade e de interesse. Em zonas muito estreitas o
  app reduz o raio e avisa.

## Como o tráfego de dados é reduzido
Busca STAC só com metadados e uma busca por fazenda; triagem de nuvem pela banda SCL (20 m) antes de baixar vermelho/NIR;
orçamento de cenas boas por janela (padrão 6); leitura de recortes COG do talhão (sem cena inteira), por overview a 20 m;
RGB de uma única cena; cache em disco por cena (reanálises não baixam nada); sem tiles de mapa externos.
Na simulação, 114 de 614 cenas tiveram NDVI lido e as leituras caíram ~81% frente a ler a série inteira.

## O que foi e o que não foi verificado
- Testado: motor, KML (relido pelo GDAL), cache, entradas, painel/mapa, app (Streamlit AppTest), leitura de COG com
  overview, fallback sem overview, e os 14 contornos reais (com dados sintéticos).
- **Não testado: a busca real no Earth Search/AWS** (sem rede no ambiente de desenvolvimento) e o aceite do KML pela
  plataforma EarthOptics. Faça primeiro um talhão pequeno e envie o erro, se houver.
- As camadas e o zoneamento seguem a hipótese de que vigor persistente reflete solo. Ela ainda precisa ser validada
  contra os resultados de TruBio (`zonas_persistencia.py` extrai as métricas por célula para esse teste).
- O Copernicus GLO-30 é modelo de superfície de 30 m: em talhões muito planos o relevo pesa pouco.

## Executar pelo GitHub (versão organizada)

Este diretório deve ser a **raiz do repositório** (não envie apenas o ZIP dentro do GitHub).

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements_zonas.txt
streamlit run app_zonas.py
```

Abra o endereço local exibido no terminal (geralmente `http://localhost:8501`).
Para validar antes de publicar:

```bash
python -m pip install pytest
python -m pytest -q tests
python zonas_core.py --selftest
```

O workflow `.github/workflows/ci.yml` executa os testes automaticamente em pushes,
Pull Requests ou manualmente em **Actions > CI · Zonas TruBio > Run workflow**.
Ele **não executa análises reais via satélite automaticamente**, não grava dados de
fazendas no repositório e não publica resultados. Não configure nenhuma credencial
no código ou em commits.

### Segurança de arquivos

O importador agora bloqueia entradas ZIP/KMZ com diretórios ascendentes, caminhos
absolutos, links simbólicos, arquivos individuais muito grandes e pacotes
excessivamente grandes. Limites padrão: 500 entradas, 150 MB por item e
500 MB descompactados por arquivo. Uploads diretos no app têm limite explícito
de 150 MB por arquivo, além dos limites configurados no servidor Streamlit.

### Implantação e limitações

- Para hospedar na Streamlit Community Cloud, aponte o arquivo principal para `app_zonas.py`.
  Dependências geoespaciais podem exigir ajustes na imagem/ambiente da hospedagem.
- O cache local pode ser efêmero no servidor; uma hospedagem persistente ou armazenamento
  externo será necessário para conservar o cache entre reinicializações.
- Os testes offline não verificam disponibilidade da API STAC, qualidade das imagens
  reais, tempo/custo de transferência nem o aceite efetivo de KML por sistemas terceiros.
- Não interprete zonas de baixo vigor automaticamente como deficiência de solo:
  clima, pragas, manejo, estágio fenológico e variedade também podem influenciar NDVI.
