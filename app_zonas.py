# -*- coding: utf-8 -*-
"""
Zonas de amostragem direcionada · TruBio  (aplicativo interno)

    pip install -r requirements_zonas.txt
    streamlit run app_zonas.py

Entrada : KML/KMZ, shapefile (.zip ou .shp+.dbf+.shx+.prj) ou GeoJSON com o limite do(s) talhão(ões).
Saída   : KML das zonas propostas para upload na plataforma EarthOptics (um polígono = uma amostra composta).
Fontes  : Sentinel-2 e Copernicus DEM via Earth Search (AWS), públicas e gratuitas. Veja zonas_core.py
          para as medidas de redução de tráfego de dados.
"""
import tempfile
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

import zonas_core as zc
from arquivos_seguros import extrair_zip_seguro


def mostrar_html(html, altura):
    """st.iframe nas versões novas do Streamlit; components.html nas antigas."""
    if hasattr(st, "iframe"):
        st.iframe(html, height=altura)
    else:
        components.html(html, height=altura)


st.set_page_config(page_title="Zonas direcionadas · TruBio", layout="wide")

EXT = {".shp", ".kml", ".kmz", ".geojson", ".json", ".gpkg"}


def coletar(arquivos_up):
    """Grava os uploads em disco (shapefile = várias partes) e devolve os caminhos legíveis."""
    d = Path(tempfile.mkdtemp(prefix="up_"))
    for f in arquivos_up:
        conteudo = f.getvalue()
        if len(conteudo) > 150 * 1024 * 1024:
            raise ValueError(f"Upload muito grande: {f.name} (máximo 150 MB por arquivo).")
        (d / Path(f.name).name).write_bytes(conteudo)
    for z in list(d.glob("*.zip")):
        extrair_zip_seguro(z, d / z.stem)
    return sorted(p for p in d.rglob("*") if p.suffix.lower() in EXT)


# ------------------------------------------------------------------ barra lateral
sb = st.sidebar
sb.header("1 · Talhões")
ups = sb.file_uploader("Limite do talhão (KML, KMZ, shapefile .zip, GeoJSON)", accept_multiple_files=True,
                       type=["kml", "kmz", "geojson", "json", "zip", "shp", "shx", "dbf", "prj", "cpg"], key="up")
unir = sb.checkbox("Os polígonos de cada arquivo são células de um grid: unir em um talhão por arquivo",
                   value=False, key="unir",
                   help="Marque se o arquivo é um grid de amostragem existente. Deixe desmarcado se cada polígono é um talhão.")
demo = sb.checkbox("Modo demonstração (dados SINTÉTICOS, sem rede)", value=False, key="demo",
                   help="Gera imagens e relevo artificiais sobre os seus polígonos (ou um pivô de exemplo). Só para conhecer a ferramenta.")

sb.header("2 · Amostragem")
dens = sb.number_input("Densidade-alvo (ha por amostra)", 4.0, 60.0, 12.0, 1.0, key="dens")
c1, c2 = sb.columns(2)
amin = c1.number_input("Área mín. (ha)", 1.0, 60.0, 6.0, 1.0, key="amin")
amax = c2.number_input("Área máx. (ha)", 2.0, 120.0, 20.0, 1.0, key="amax")
rig = sb.slider("Rigidez dos limites de área", 0.0, 0.3, 0.05, 0.01, key="rig",
                help="Os limites são FLEXÍVEIS: 0 ignora; valores maiores forçam mais, ao custo de zonas menos homogêneas.")
comp = sb.slider("Compacidade das zonas", 0.0, 1.0, 0.5, 0.1, key="comp",
                 help="0 = zonas seguem só o padrão de vigor/relevo; 1 = formas mais compactas e fáceis de percorrer.")
suav = sb.slider("Suavização dos contornos", 0.0, 1.0, 0.8, 0.1, key="suav",
                 help="0 = contornos em degraus de pixel; valores maiores arredondam zonas e camadas (as fronteiras "
                      "compartilhadas são suavizadas uma só vez, sem criar frestas nem sobreposições).")

with sb.expander("Dados de satélite e tráfego"):
    res_m = st.selectbox("Resolução de trabalho (m)", [20, 10], index=0, key="resolucao",
                         help="20 m lê o nível reduzido do COG: menos dados e mais rápido. 10 m detalha manchas pequenas.")
    obs = st.number_input("Cenas boas por janela de safra", 3, 20, 6, 1, key="obs",
                          help="Orçamento de download: para de baixar quando já há esta quantidade de cenas sem nuvem.")
    nuv = st.slider("Nuvem máxima da cena (metadados)", 10, 100, 60, 5, key="nuv")
    claro = st.slider("Fração do talhão sem nuvem para aceitar a cena", 0.3, 1.0, 0.6, 0.05, key="claro")
    ini = st.text_input("Início da série", "2017-01-01", key="ini")
    buf = st.number_input("Recuo da borda (m)", 0.0, 100.0, 20.0, 5.0, key="buf")
    cache_dir = st.text_input("Pasta de cache", "cache_zonas", key="cache")

with sb.expander("Camadas: estabilidade e interesse"):
    lim = st.slider("Desvio relativo à mediana considerado relevante", 0.05, 0.30, 0.10, 0.01, key="lim")
    f_int = st.slider("Queda persistente: janelas abaixo do limiar (≥)", 0.2, 0.9, 0.5, 0.05, key="fint")
    f_est = st.slider("Estável: janelas com desvio (≤)", 0.0, 0.4, 0.15, 0.05, key="fest")
    nrec = st.number_input("Janelas recentes avaliadas", 2, 6, 3, 1, key="nrec")
    min_cam = st.number_input("Área mínima de uma camada (ha)", 0.2, 10.0, 1.0, 0.2, key="mincam")
    peso_rel = st.slider("Peso do relevo no zoneamento", 0.0, 1.5, 0.5, 0.1, key="prel")

rodar = sb.button("▶ Analisar", type="primary", key="rodar")

# ------------------------------------------------------------------ cabeçalho
st.title("Zonas de amostragem direcionada")
st.caption("Grid variável a partir do que a planta mostra ao longo dos anos (NDVI Sentinel-2) e do relevo (Copernicus DEM). "
           "Cada zona = uma amostra composta. Fontes públicas e gratuitas, com tráfego de dados minimizado.")

if rodar:
    if amin > amax:
        st.error("A área mínima não pode superar a área máxima.")
        st.stop()
    try:
        from datetime import date
        data_inicio = date.fromisoformat(ini)
        if data_inicio > date.today():
            raise ValueError("A data inicial não pode estar no futuro.")
    except ValueError as e:
        st.error(f"Data inicial inválida (use AAAA-MM-DD): {e}")
        st.stop()
    p = zc.Params(densidade_ha=dens, area_min_ha=amin, area_max_ha=amax, rigidez_area=rig, compacidade=comp, suavizacao=suav,
                  res_m=int(res_m), obs_por_janela=int(obs), nuvem_max_cena=float(nuv), claro_min=claro,
                  inicio=ini, buffer_borda_m=buf, cache_dir=cache_dir, limiar_anomalia=lim, f_interesse=f_int,
                  f_estavel=f_est, n_recente=int(nrec), area_min_camada_ha=min_cam, peso_relevo=peso_rel)
    try:
        if ups:
            talhoes = zc.ler_talhoes(coletar(ups), unir_por_arquivo=unir)
        elif demo:
            talhoes = zc.talhao_exemplo()
        else:
            talhoes = []
    except Exception as e:
        talhoes = []
        st.error(f"Não consegui ler o(s) arquivo(s): {e}")
    if not talhoes:
        st.warning("Envie o limite de um talhão (ou marque o modo demonstração).")
    else:
        fonte = zc.FonteSintetica() if demo else zc.FonteSTAC()
        barra, msg = st.progress(0.0), st.empty()
        erros = []

        def prog(f, m):
            barra.progress(float(min(max(f, 0.0), 1.0)))
            msg.write(m)
        res = zc.analisar_fazenda(talhoes, fonte, p, prog, erros)
        barra.empty(); msg.empty()
        st.session_state["_resultados"], st.session_state["_erros"] = res, erros

res = st.session_state.get("_resultados")
if not res:
    st.info("Envie o limite do talhão na barra lateral, ajuste a densidade e clique em **Analisar**. "
            "Use o **modo demonstração** para conhecer a ferramenta sem baixar imagens.")
    with st.expander("Como funciona"):
        st.markdown(
            "- **Série:** para cada safra (verão: out–mar; inverno/safrinha: abr–set) calcula-se o NDVI de pico por pixel e o "
            "desvio relativo à **mediana do talhão**. Janelas sem lavoura são descartadas.\n"
            "- **Zoneamento:** pixels contíguos com histórico e relevo parecidos são agrupados em N = área ÷ densidade zonas, "
            "com limites de área flexíveis e formas compactas.\n"
            "- **Estabilidade:** áreas sem desvios relevantes sobre a mediana em (quase) todas as janelas.\n"
            "- **Interesse:** queda persistente (abaixo da mediana na maioria das janelas) ou queda recente "
            "(abaixo nas últimas janelas e normal antes).\n"
            "- **Tráfego:** triagem de nuvem pela banda SCL antes de baixar vermelho/NIR, orçamento de cenas por janela, "
            "recortes COG do talhão e cache em disco.")
    st.stop()

for nome, e in st.session_state.get("_erros", []):
    st.error(f"**{nome}** não foi analisado: {e}")

if res[0].fonte.startswith("SINT"):
    st.warning("DADOS SINTÉTICOS: imagens e relevo artificiais, só para demonstração. Não use estes resultados para pedir análises.")

# ------------------------------------------------------------------ downloads
st.subheader("Saída para a plataforma EarthOptics")
d1, d2 = st.columns([2, 3])
limpo = d1.checkbox("KML limpo (só nome e geometria)", True, key="limpo",
                    help="Recomendado para upload. Desmarque para incluir cores e descrição (bom para Google Earth).")
padrao = d2.text_input("Padrão do nome da zona", "{talhao}_Z{n:02d}", key="padrao",
                       help="Variáveis: {talhao} e {n}. Ex.: {talhao}_Z{n:02d}")
try:
    kml = zc.kml_zonas(res, limpo=limpo, padrao_nome=padrao)
except Exception as e:
    st.error(f"Padrão de nome inválido: {e}")
    kml = zc.kml_zonas(res, limpo=limpo)
b = st.columns(4)
b[0].download_button("⬇ KML para upload (zonas)", kml.encode("utf-8"), "zonas_propostas.kml",
                     "application/vnd.google-earth.kml+xml", key="dl_kml")
b[1].download_button("⬇ KML de camadas (estabilidade/interesse)", zc.kml_camadas(res).encode("utf-8"),
                     "camadas_analise.kml", "application/vnd.google-earth.kml+xml", key="dl_cam")
b[2].download_button("⬇ Tabela de zonas (CSV)", zc.tabela_zonas(res).to_csv(index=False).encode("utf-8"),
                     "zonas.csv", "text/csv", key="dl_csv")

# ------------------------------------------------------------------ visão por talhão
st.divider()
nomes = [r.nome for r in res]
sel = st.selectbox("Talhão", nomes, key="sel") if len(nomes) > 1 else nomes[0]
r = res[nomes.index(sel)]
s = r.resumo()
m = st.columns(5)
m[0].metric("Zonas (amostras)", s["n_zonas"])
m[1].metric("ha por amostra", f"{s['ha_por_amostra']:.1f}")
m[2].metric("Área das zonas", f"{s['area_min']:.1f} – {s['area_max']:.1f} ha", f"{s['fora_dos_limites']} fora dos limites", delta_color="off")
m[3].metric("Estabilidade", f"{s['pct_estavel']:.0f}% da área")
m[4].metric("Interesse (queda de vigor)", f"{s['pct_interesse']:.0f}% da área")
for a in r.avisos:
    st.caption(f"⚠ {a}")

t_mapa, t_painel, t_zonas, t_dados = st.tabs(["Mapa interativo", "Painel (4 vistas)", "Zonas", "Dados baixados"])
with t_mapa:
    st.caption("Use o controle de camadas (canto superior direito): RGB, clusterização, zonas de estabilidade e zonas de interesse.")
    mostrar_html(zc.mapa_html(r), 700)
with t_painel:
    png = zc.painel_png(r)
    st.image(png)
    st.download_button("⬇ Painel (PNG)", png, f"painel_{zc._slug(r.nome)}.png", "image/png", key="dl_png")
with t_zonas:
    t = r.tabela()
    mostrar = t[["nome", "area_ha", "pct_estavel", "pct_interesse", "mu", "fneg", "elev_rel", "decliv_pct", "twi"]].rename(columns={
        "nome": "Zona", "area_ha": "Área (ha)", "pct_estavel": "% estável", "pct_interesse": "% interesse", "mu": "Desvio médio",
        "fneg": "Fração de janelas em queda", "elev_rel": "Elev. relativa (m)", "decliv_pct": "Declividade (%)", "twi": "TWI"})
    st.dataframe(mostrar.round(2), hide_index=True)
    st.bar_chart(t.set_index("nome")["area_ha"])
with t_dados:
    tr = r.trafego
    ingenuo = tr.cenas_encontradas * 3          # ler vermelho + NIR + SCL de TODAS as cenas
    nec = tr.leituras_necessarias
    st.markdown(
        f"**Fonte:** {r.fonte} · **Janelas de safra usadas:** {len(r.janelas_usadas)} de {len(r.janelas_total)} "
        f"({', '.join(r.janelas_usadas)})")
    c = st.columns(4)
    c[0].metric("Cenas encontradas (só metadados)", tr.cenas_encontradas)
    c[1].metric("Triadas pela SCL (rede / cache)", f"{tr.triagens_rede} / {tr.triagens_cache}")
    c[2].metric("NDVI baixado (rede / cache)", f"{tr.ndvi_rede} / {tr.ndvi_cache}")
    c[3].metric("Descartadas por nuvem", tr.cenas_descartadas_nuvem)
    st.markdown(
        f"**Leituras de recorte que esta configuração exige:** {nec} (cena × banda), contra cerca de {ingenuo} se lesse "
        f"vermelho + NIR + SCL de toda a série: **{100 * (1 - nec / max(ingenuo, 1)):.0f}% a menos**. "
        f"Nesta execução foram {tr.leituras_banda} pela rede; o restante já estava no cache.")
    if tr.bytes_medidos is not None:
        st.caption(f"Tráfego de rede medido durante a aquisição (aproximado, inclui outros programas): {tr.bytes_medidos / 1e6:.1f} MB.")
    else:
        st.caption("Instale `psutil` para medir o tráfego de rede real durante a análise.")
    st.caption(f"RGB: cena de {r.rgb_data or '—'} (uma única cena, a mais limpa dos últimos {r.params.dias_rgb} dias).")
