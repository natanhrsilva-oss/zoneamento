# Relatório de revisão — TruBio / Zonas de amostragem
**Data:** 30/10/2026

## 1. Material analisado
O pacote original tinha quatro arquivos: duas prévias PNG, um HTML demonstrativo e
`zonas_app.zip`. O aplicativo real estava dentro do ZIP aninhado, com os módulos
`app_zonas.py`, `zonas_core.py`, `zonas_persistencia.py`, `zonas_selftest.py`, o
arquivo de dependências e a documentação. As prévias e o HTML foram preservados
no pacote de entrega, em `exemplos/`.

## 2. Melhorias aplicadas

1. **Segurança de arquivos compactados:** criei `arquivos_seguros.py`, que substitui
   `ZipFile.extractall` nos três pontos de importação (app, motor e módulo de
   persistência). Bloqueia `../`, caminhos absolutos, letra de unidade Windows,
   links simbólicos e limites excessivos de tamanho/quantidade de entradas.
   A leitura ocorre em blocos e valida o total efetivamente extraído.
2. **Testes de regressão de segurança:** `tests/test_arquivos_seguros.py` cobre
   ZIP válido, quatro caminhos perigosos e tamanho superior ao permitido.
3. **CI no GitHub:** `.github/workflows/ci.yml` usa Python 3.11, instala dependências,
   compila os módulos e executa os testes adicionais e o selftest agronômico offline.
   Há timeout de 25 minutos, concorrência controlada e permissão apenas de leitura.
4. **Experiência de uso:** o formulário verifica se a área mínima não supera a
   máxima e se a data inicial está em ISO AAAA-MM-DD e não está no futuro.
5. **Higiene do repositório:** `.gitignore` exclui caches, ambientes locais,
   segredos e arquivos transitórios; `.python-version` declara Python 3.11.
6. **Documentação:** adicionei instruções de instalação, teste e GitHub Actions,
   cuidados de hospedagem e ressalvas para interpretação do NDVI.

## 3. O que foi preservado
Mantive o motor de zoneamento, os pesos, a lógica de NDVI, a avaliação de relevo,
a suavização das fronteiras, o cache, as saídas KML/CSV/PNG e o modo demonstração.
Não alterei fórmulas agronômicas nem faixas de classificação que precisariam de
validação de campo.

## 4. Testes realizados neste ambiente
- `python -m compileall -q .`: **aprovado**.
- `python -m pytest -q tests`: **6/6 aprovados**.
- `python zonas_core.py --selftest`: **aprovado**, mensagem final
  `TODOS OS TESTES PASSARAM`.
- Rede STAC/AWS em tempo real: **não testada**.
- Interface Streamlit em navegador: **não testada** neste ambiente (pacote
  Streamlit não disponível localmente).
- Execução do workflow dentro do GitHub: **não executada**, apenas configurada.
- Aceite do KML na plataforma EarthOptics: **não testado**.

## 5. Melhorias futuras recomendadas
**Prioridade alta:** validar uma área real em diversos estágios da safra e testar o
KML na plataforma de destino; proteger a ferramenta caso seja hospedada publicamente
(com autenticação e controles de consumo de recursos); adicionar logs estruturados.

**Prioridade média:** fixar versões das bibliotecas com arquivo de lock gerado em
ambiente de CI; testar instalação limpa em Linux e Windows; criar testes de
integração da busca STAC com mocks e tolerância a falhas de serviço; monitorar o
espaço ocupado pelo cache e implementar prazo de retenção.

**Prioridade agronômica:** confrontar zonas propostas com resultados de análise de
solo/TruBio, verificar correlação e taxas de falsos positivos por cultura, período
fenológico, manejo, irrigação e índice de vegetação; só então recalibrar thresholds.

## 6. Publicação
Crie um repositório no GitHub e envie o **conteúdo da pasta `zonas_app/` para a
raiz**, incluindo `.github/workflows/ci.yml`. Confira a aba **Actions** depois do
primeiro push. Execute `streamlit run app_zonas.py` na máquina ou configure uma
hospedagem Python adequada. Nunca publique dados privados, `.env` ou arquivos de
credenciais.

**Resumo:** versão preparada para GitHub com melhorias de segurança e testes
locais aprovados, mas não homologada ainda contra fontes externas reais.
