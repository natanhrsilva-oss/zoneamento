# -*- coding: utf-8 -*-
"""Testes automáticos (sem rede) do motor de zoneamento. Uso: python zonas_core.py --selftest"""
import io
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

import zonas_core as zc
import zonas_persistencia as zp
from shapely.geometry import Point, Polygon

_ok = True


def check(nome, cond, extra=""):
    global _ok
    _ok &= bool(cond)
    print(f"[{'OK' if cond else 'FALHA'}] {nome} {extra}")


def pivo(lon=-47.9, lat=-14.7, raio_m=440.0):
    """Círculo de ~60 ha em WGS84 (construído em plano local, sem UTM)."""
    c = zp.crs_local(Point(lon, lat))
    g = Point(0, 0).buffer(raio_m, 64)
    return zp.shp_transform(zp.func_transf(c, "EPSG:4326"), g)


def t_zoneamento():
    rng = np.random.default_rng(0)
    ny = nx = 60
    verdade = np.zeros((ny, nx), int)
    verdade[:, 20:40] = 1
    verdade[:, 40:] = 2
    f = np.stack([verdade * 2.0, (verdade == 1) * 1.5, np.zeros((ny, nx))], -1) + rng.normal(0, 0.35, (ny, nx, 3))
    mask = np.ones((ny, nx), bool)
    p = zc.Params(densidade_ha=48.0, area_min_ha=10, area_max_ha=100, piso_ha=2.0)
    px_ha = 0.04
    ss = float(np.sum((f.reshape(-1, 3) - f.reshape(-1, 3).mean(0)) ** 2))
    lab = zc.zonear(f, zc.xy_padronizado(mask), mask, px_ha, p, ss)
    n = len(np.unique(lab))
    # pureza: cada zona deve cair majoritariamente numa região verdadeira
    pur = sum(np.bincount(verdade[lab == k]).max() for k in np.unique(lab)) / lab.size
    check("zoneamento recupera 3 regiões (3 zonas, pureza > 0,95)", n == 3 and pur > 0.95, f"n={n} pureza={pur:.3f}")


def t_limites_area():
    rng = np.random.default_rng(1)
    ny = nx = 100
    f = np.stack([ndi_suave(rng, ny, nx) for _ in range(4)], -1)
    mask = np.ones((ny, nx), bool)
    ss = float(np.sum((f.reshape(-1, 4) - f.reshape(-1, 4).mean(0)) ** 2))
    res = {}
    for rig in (0.0, 0.2):
        p = zc.Params(densidade_ha=20, area_min_ha=12, area_max_ha=30, rigidez_area=rig, piso_ha=2)
        lab = zc.zonear(f, zc.xy_padronizado(mask), mask, 0.04, p, ss)
        a = np.bincount(lab.ravel())[1:] * 0.04
        res[rig] = (len(a), float(np.mean((a >= 12) & (a <= 30))), a.min(), a.max())
    check("limites de área flexíveis: rigidez 0,2 coloca mais zonas dentro de [12, 30] ha que rigidez 0",
          res[0.2][1] > res[0.0][1] and res[0.2][1] >= 0.7,
          f"livre={res[0.0][1]:.0%} (min {res[0.0][2]:.1f}, max {res[0.0][3]:.1f}) | rígido={res[0.2][1]:.0%} (min {res[0.2][2]:.1f}, max {res[0.2][3]:.1f})")
    check("número de zonas igual ao alvo (400 ha / 20 = 20)", res[0.0][0] == 20 and res[0.2][0] == 20, f"{res[0.0][0]}, {res[0.2][0]}")


def ndi_suave(rng, ny, nx):
    from scipy import ndimage as ndi
    a = ndi.gaussian_filter(rng.normal(size=(ny, nx)), 9)
    return a / a.std()


def t_pipeline():
    g = pivo(raio_m=620)
    area_campo = zp.area_geodesica_ha(g)
    cache = tempfile.mkdtemp(prefix="zc_cache_")
    p = zc.Params(cache_dir=cache, inicio="2017-01-01", fim="2026-10-30")
    r = zc.analisar_talhao("Pivo_teste", g, zc.FonteSintetica(), p)
    z = r.zonas
    soma = sum(x["area_ha"] for x in z)
    check("zonas cobrem o talhão (±1,5%)", abs(soma - area_campo) / area_campo < 0.015, f"{soma:.2f} de {area_campo:.2f} ha")
    uni = zp.unary_union([x["geom_ll"] for x in z])
    check("sem sobreposição entre zonas", abs(sum(x["geom_ll"].area for x in z) - uni.area) / uni.area < 0.002)
    check("polígonos válidos", all(x["geom_ll"].is_valid for x in z))
    alvo = round(area_campo / p.densidade_ha)
    check(f"nº de zonas igual ao alvo ({alvo}) ±1", abs(len(z) - alvo) <= 1, f"{len(z)} zonas")
    a = np.array([x["area_ha"] for x in z])
    dentro = np.mean((a >= p.area_min_ha) & (a <= p.area_max_ha))
    check("limites flexíveis: >= 60% das zonas dentro de [mín, máx] e nenhuma além de 50% dos limites",
          dentro >= 0.6 and a.min() >= 0.5 * p.area_min_ha and a.max() <= 1.5 * p.area_max_ha,
          f"{dentro:.0%} dentro | {[round(x, 1) for x in a]}")
    t = r.trafego
    jan = len(r.janelas_total)
    check("orçamento de cenas respeitado (<= janelas x obs_por_janela)", t.ndvi_rede <= jan * p.obs_por_janela, f"{t.ndvi_rede} <= {jan * p.obs_por_janela}")
    total_cenas = t.cenas_encontradas
    check("baixa uma fração pequena da série", t.ndvi_rede < 0.35 * total_cenas, f"{t.ndvi_rede} NDVI de {total_cenas} cenas encontradas")
    check("janelas de pousio descartadas (usadas < total)", len(r.janelas_usadas) < len(r.janelas_total), f"{len(r.janelas_usadas)}/{len(r.janelas_total)}")
    check("camada de estabilidade e de interesse não vazias", len(r.estavel) > 0 and len(r.interesse) > 0,
          f"estável={len(r.estavel)} interesse={len(r.interesse)}")
    tipos = {x["tipo"] for x in r.interesse}
    check("detecta queda persistente e queda recente", {"queda persistente", "queda recente"} <= tipos, f"{tipos}")
    # as quedas sintéticas ficam a oeste (persistente) e a leste/sul (recente): conferir lado
    cx = g.centroid.x
    pers = [x for x in r.interesse if x["tipo"] == "queda persistente"]
    rec = [x for x in r.interesse if x["tipo"] == "queda recente"]
    check("queda persistente a oeste e recente a leste do centro (posição da mancha sintética)",
          np.mean([x["geom_ll"].centroid.x for x in pers]) < cx < np.mean([x["geom_ll"].centroid.x for x in rec]))
    ue, ui = zp.unary_union([x["geom_ll"] for x in r.estavel]), zp.unary_union([x["geom_ll"] for x in r.interesse])
    check("camadas de estabilidade e de interesse não se sobrepõem", ue.intersection(ui).area < 1e-4 * ui.area, f"sobreposição={ue.intersection(ui).area / ui.area:.2e}")
    check("RGB gerada (3 bandas, uint8)", r.rgb is not None and r.rgb.shape[0] == 3 and r.rgb.dtype == np.uint8)
    # segunda execução: tudo vem do cache
    r2 = zc.analisar_talhao("Pivo_teste", g, zc.FonteSintetica(), p)
    t2 = r2.trafego
    check("reanálise sem baixar nada (tudo do cache)", t2.ndvi_rede == 0 and t2.triagens_rede == 0 and t2.rgb_rede == 0,
          f"ndvi_rede={t2.ndvi_rede} triagens_rede={t2.triagens_rede} rgb_rede={t2.rgb_rede} | ndvi_cache={t2.ndvi_cache}")
    check("reanálise reproduz o mesmo zoneamento", len(r2.zonas) == len(z) and all(abs(a["area_ha"] - b["area_ha"]) < 1e-6 for a, b in zip(z, r2.zonas)))
    return r


def t_kml(r):
    k = zc.kml_zonas([r])
    raiz = ET.fromstring(k)
    ns = {"k": "http://www.opengis.net/kml/2.2"}
    pm = raiz.findall(".//k:Placemark", ns)
    check("KML bem-formado com 1 Placemark por zona", len(pm) == len(r.zonas), f"{len(pm)} placemarks")
    check("KML limpo: só nome e polígono", all(len(list(x)) == 2 for x in pm))
    nomes = [x.find("k:name", ns).text for x in pm]
    check("nomes únicos", len(set(nomes)) == len(nomes), f"{nomes[:3]}...")
    tmp = Path(tempfile.mkdtemp()) / "z.kml"
    tmp.write_text(k, encoding="utf-8")
    import geopandas as gpd
    g = gpd.read_file(tmp, driver="KML")
    a_rel = sum(zp.area_geodesica_ha(x) for x in g.geometry)
    a_orig = sum(z["area_ha"] for z in r.zonas)
    check("KML relido pelo GDAL: mesma contagem e área", len(g) == len(r.zonas) and abs(a_rel - a_orig) / a_orig < 1e-4, f"{len(g)} zonas, {a_rel:.2f} ha")
    def ccw(txt):
        pts = [tuple(map(float, q.split(",")[:2])) for q in txt.split()]
        return sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in zip(pts, pts[1:])) > 0
    ext = [x.text for x in raiz.findall(".//k:outerBoundaryIs/k:LinearRing/k:coordinates", ns)]
    check("anel externo em sentido anti-horário no KML (padrão OGC)", len(ext) >= len(r.zonas) and all(ccw(c) for c in ext), f"{len(ext)} anéis")
    kc = zc.kml_camadas([r])
    ET.fromstring(kc)
    check("KML de camadas bem-formado", "Zonas de estabilidade" in kc and "queda persistente" in kc.lower())
    kf = zc.kml_zonas([r], limpo=False)
    ET.fromstring(kf)
    check("KML completo com estilos bem-formado", "<Style" in kf and "<description>" in kf)
    return tmp


def t_entrada(kml_path, r):
    import geopandas as gpd, zipfile
    d = Path(tempfile.mkdtemp())
    gj = d / "talhao.geojson"
    gpd.GeoDataFrame({"Name": [r.nome]}, geometry=[r.campo_ll], crs=4326).to_file(gj, driver="GeoJSON")
    shp = d / "Talhao_SHP.shp"
    gpd.GeoDataFrame({"talhao": ["A1"]}, geometry=[r.campo_ll], crs=4326).to_file(shp)
    kmz = d / "t.kmz"
    with zipfile.ZipFile(kmz, "w") as z:
        z.write(kml_path, "doc.kml")
    a = zc.ler_talhoes([gj, shp, kmz])
    check("lê GeoJSON, shapefile e KMZ", len(a) >= 3, f"{[n for n, _ in a][:6]}")
    check("nomes vêm dos atributos (Name / talhao)", a[0][0] == r.nome and a[1][0] == "A1")
    u = zc.ler_talhoes([kml_path], unir_por_arquivo=True)
    check("'unir por arquivo' dissolve as zonas num talhão só", len(u) == 1 and abs(zp.area_geodesica_ha(u[0][1]) - zp.area_geodesica_ha(r.campo_ll)) / zp.area_geodesica_ha(r.campo_ll) < 0.015)


def t_leitura_cog():
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_origin
    n, x0, y0 = 2048, 600000.0, 8850000.0
    base = (np.random.default_rng(3).random((n, n)) * 3000 + 2000).astype("uint16")
    d = Path(tempfile.mkdtemp())
    prof = dict(driver="GTiff", height=n, width=n, count=1, dtype="uint16", crs="EPSG:31981", transform=from_origin(x0, y0, 10, 10),
                tiled=True, blockxsize=512, blockysize=512, compress="deflate", nodata=0)
    with rasterio.open(d / "cog.tif", "w", **prof) as f:
        f.write(base, 1); f.build_overviews([2, 4], Resampling.average)
    with rasterio.open(d / "plain.tif", "w", **{**prof, "tiled": False}) as f:
        f.write(base, 1)
    gr = zc.Grade(31981, 600400.0, 8850000.0 - 1000, 20.0, 40, 40)
    ref = base[100 + 0:100 + 80, 40:120].astype(float).reshape(40, 2, 40, 2).mean((1, 3))
    a = zc.ler_janela(str(d / "cog.tif"), gr, Resampling.average, 0)
    b = zc.ler_janela(str(d / "plain.tif"), gr, Resampling.average, 0)      # sem overview: cai para decimação
    check("leitura via overview do COG = média 2x2 (±1 DN)", np.abs(a - ref).max() <= 1.0, f"max|dif|={np.abs(a - ref).max():.1f}")
    check("arquivo sem overview: fallback funciona", np.abs(b - ref).max() <= 1.0)
    fora = zc.Grade(31981, 599600.0, 8850000.0, 20.0, 40, 40)   # cruza a borda oeste da cena (x=600000)
    c = zc.ler_janela(str(d / "cog.tif"), fora, Resampling.average, 0)
    check("janela parcialmente fora da cena: preenche com 0 (nodata)", (c == 0).any() and (c > 0).any())


def t_visual(r):
    png = zc.painel_png(r, dpi=60)
    check("painel PNG com as 4 vistas", png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 20000, f"{len(png) / 1024:.0f} kB")
    h = zc.mapa_html(r)
    check("mapa interativo: camadas RGB/zonas/estabilidade/interesse", all(s in h for s in ("Clusteriza", "Zonas de estabilidade", "Zonas de interesse", "RGB (")))
    check("mapa sem tiles externos de mapa-base", "tile.openstreetmap" not in h and "arcgisonline" not in h and "cartocdn" not in h)


def t_suavizacao():
    """A suavização não pode criar frestas/sobreposições e deve reduzir os contornos em degraus."""
    g = pivo(raio_m=620)
    area_campo = zp.area_geodesica_ha(g)
    cache = tempfile.mkdtemp(prefix="zc_cache_")
    r0 = zc.analisar_talhao("liso0", g, zc.FonteSintetica(), zc.Params(cache_dir=cache, suavizacao=0.0))
    r1 = zc.analisar_talhao("liso1", g, zc.FonteSintetica(), zc.Params(cache_dir=cache, suavizacao=0.8))

    def eixo(r):                       # fração do perímetro das zonas em segmentos horizontais/verticais (degraus)
        tot = al = 0.0
        for z in r.zonas:
            for q in zc._poligonos_de(z["geom_utm"]):
                d = np.diff(np.asarray(q.exterior.coords), axis=0)
                L = np.hypot(d[:, 0], d[:, 1])
                tot += L.sum(); al += L[(np.abs(d[:, 0]) < 1e-6) | (np.abs(d[:, 1]) < 1e-6)].sum()
        return al / tot
    e0, e1 = eixo(r0), eixo(r1)
    check("suavização reduz os contornos em degraus (perímetro alinhado aos eixos)", e1 < 0.4 * e0, f"{e0:.0%} -> {e1:.0%}")
    soma = sum(z["area_ha"] for z in r1.zonas)
    u = zp.unary_union([z["geom_ll"] for z in r1.zonas])
    check("suavização mantém cobertura de 100% (±0,5%)", abs(soma - area_campo) / area_campo < 0.005, f"{soma / area_campo:.2%}")
    check("suavização não cria sobreposição entre zonas", abs(sum(z["geom_ll"].area for z in r1.zonas) - u.area) / u.area < 1e-5)
    check("suavização mantém geometrias válidas e o mesmo nº de zonas",
          all(z["geom_ll"].is_valid for z in r1.zonas) and len(r1.zonas) == len(r0.zonas), f"{len(r0.zonas)} -> {len(r1.zonas)}")
    a0, a1 = sorted(z["area_ha"] for z in r0.zonas), sorted(z["area_ha"] for z in r1.zonas)
    check("áreas das zonas quase inalteradas (mín e máx ±10%)", abs(a1[0] - a0[0]) / a0[0] < 0.10 and abs(a1[-1] - a0[-1]) / a0[-1] < 0.10,
          f"mín {a0[0]:.1f}->{a1[0]:.1f}, máx {a0[-1]:.1f}->{a1[-1]:.1f}")
    ue, ui = zp.unary_union([x["geom_ll"] for x in r1.estavel]), zp.unary_union([x["geom_ll"] for x in r1.interesse])
    check("camadas suavizadas continuam sem sobreposição", ue.intersection(ui).area < 1e-4 * ui.area)


def t_rede_lacos():
    """Regressão: zona-ilha e laço preso a um nó compartilhado não podem desprender a rede."""
    from shapely.geometry import Polygon as P, box
    A = box(0, 0, 1000, 2000)
    S = P([(1000, 1000), (1300, 900), (1600, 1000), (1300, 1100)])         # losango (6 ha) com um vértice no nó (1000, 1000)
    B = box(1000, 0, 2000, 2000).difference(S)                              # tem um furo encostado num ponto
    I = box(1500, 1500, 1800, 1800)                                         # ilha isolada (9 ha) dentro de B
    B = B.difference(I)
    zonas = [A, B, S, I]
    uni = zp.unary_union(zonas)
    novas = zc.suavizar_rede(zonas, uni, 30.0)
    check("rede com laço preso a nó e ilha: suavização aplicada (não devolve o original)", novas is not zonas)
    ratio = [n.area / o.area for n, o in zip(novas, zonas)]
    check("áreas de cada zona preservadas (±15%)", all(0.85 < r < 1.15 for r in ratio), f"{[round(r, 3) for r in ratio]}")
    check("cobertura total preservada", abs(sum(n.area for n in novas) - uni.area) / uni.area < 0.002)
    check("geometrias válidas e sem sobreposição", all(n.is_valid for n in novas) and
          abs(sum(n.area for n in novas) - zp.unary_union(novas).area) < 1e-3 * uni.area)
    check("a borda externa do conjunto não se moveu", abs(zp.unary_union(novas).symmetric_difference(uni).area) < 0.002 * uni.area)


def t_triagem():
    """Cenas sob nuvem não devem ter NDVI baixado."""
    g = pivo(-54.85, -10.39, 300)
    cache = tempfile.mkdtemp(prefix="zc_cache_")

    class FonteNublada(zc.FonteSintetica):
        def claro(self, cena, grade20, mask20):
            return 0.1 if cena.id.endswith(("01", "02", "03", "04", "05", "06", "07", "08", "09")) else 0.95
    p = zc.Params(cache_dir=cache, obs_por_janela=4)
    r = zc.analisar_talhao("nublado", g, FonteNublada(), p)
    check("cenas com nuvem sobre o talhão são descartadas sem baixar NDVI", r.trafego.cenas_descartadas_nuvem > 0 and
          r.trafego.ndvi_rede <= len(r.janelas_total) * 4, f"descartadas={r.trafego.cenas_descartadas_nuvem}, NDVI baixados={r.trafego.ndvi_rede}")


def rodar():
    print("== Zoneamento"); t_zoneamento(); t_limites_area()
    print("== Pipeline completo (fonte sintética)"); r = t_pipeline()
    print("== KML e entradas"); k = t_kml(r); t_entrada(k, r)
    print("== Leitura de COG"); t_leitura_cog()
    print("== Visualização"); t_visual(r)
    print("== Suavização dos contornos"); t_suavizacao()
    print("== Rede de zonas (laços e ilhas)"); t_rede_lacos()
    print("== Triagem de nuvem"); t_triagem()
    print("\nTODOS OS TESTES PASSARAM" if _ok else "\nHÁ TESTES COM FALHA")
    return 0 if _ok else 1


if __name__ == "__main__":
    sys.exit(rodar())
