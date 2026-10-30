#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zonas_persistencia.py
=====================
Para cada célula (polígono) de cada talhão, extrai métricas de PERSISTÊNCIA DE NDVI
(Sentinel-2) e de RELEVO (Copernicus DEM 30 m), usando só fontes públicas e gratuitas
via STAC (Earth Search / AWS). Saída para validar se o satélite e o relevo predizem
os indicadores biológicos do TruBio dentro de cada talhão.

CONVENÇÕES
  * Entrada e saída em WGS84 (EPSG:4326). Nunca se calcula área/distância em graus:
      - área, perímetro: geodésicos (pyproj.Geod);
      - buffer de borda: plano métrico local (azimutal equidistante centrado no talhão),
        sem zona UTM a escolher;
      - imagens: lidas no CRS nativo do tile (sem reamostrar para graus); o polígono é
        levado ao CRS do tile só para ler os pixels.
  * Cada arquivo de entrada = um talhão; cada polígono dele = uma célula do grid.
    Chave de junção com o TruBio: (campo, sample_id).

MÉTODO (acordado)
  * Série 2017–hoje, duas janelas por ano: V = outubro–março (safra de verão, rotulada
    pelo ano de início) e I = abril–setembro. Só janelas completas.
  * Por pixel e janela: NDVI no percentil 90 (após máscara SCL 4/5), exigindo >= obs_min
    observações válidas.
  * Desvio relativo à média do talhão em cada janela: d = (V - média_talhão) / média_talhão.
    Janelas sem lavoura (média do talhão < --ndvi-min-janela) ou com pouca cobertura
    são descartadas.
  * Por pixel: mu (média dos d), sigma (desvio-padrão dos d), fneg / fpos (fração das
    janelas com d < -limiar / d > +limiar).
  * Relevo (Copernicus GLO-30): declividade, TPI 200 m, elevação relativa, TWI.
  * Agregação por célula (média) + variação interna da célula (_sd_intra).

INSTALAÇÃO
  pip install numpy pandas geopandas shapely pyproj scipy rasterio xarray \
              pystac-client odc-stac pyogrio

USO
  python zonas_persistencia.py PASTA_OU_ZIP_OU_SHP [...] --saida celulas_features.csv
  python zonas_persistencia.py ... --so-campo "Talhão 8"      # teste rápido num talhão
  python zonas_persistencia.py ... --so-geometria             # só valida entradas (sem rede)
  python zonas_persistencia.py --selftest                     # testes numéricos (sem rede)

LIMITES CONHECIDOS
  * O Copernicus GLO-30 é um modelo de SUPERFÍCIE (inclui dossel/construções) com
    resolução de 30 m: em talhões muito planos o sinal de relevo é fraco.
  * Células com poucas janelas válidas (n_janelas < 5) devem ser tratadas com cautela.
  * Testado com itens STAC sintéticos locais; a busca real no Earth Search depende de
    rede e deve ser conferida no primeiro talhão (use --so-campo).
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import math
import sys
import tempfile
import time
import warnings
import zipfile
from collections import Counter
from datetime import date, datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import geopandas as gpd
    import shapely
    from pyproj import CRS, Geod, Transformer
    from scipy import ndimage as ndi
    from shapely.ops import transform as shp_transform
    from shapely.ops import unary_union
    import rasterio.features
except ImportError as e:  # pragma: no cover
    sys.exit(f"Dependência ausente: {e}. Veja a seção INSTALAÇÃO no topo do arquivo.")

from arquivos_seguros import extrair_zip_seguro

STAC_URL = "https://earth-search.aws.element84.com/v1"
GEOD = Geod(ellps="WGS84")
SCL_VALIDO = (4, 5)          # vegetação e solo não vegetado
RES_S2 = 10                  # m
RES_DEM = 30                 # m
EXTENSOES = {".shp", ".kml", ".geojson", ".gpkg"}
_AVISOS: set[str] = set()


def log(*a):
    print(*a, flush=True)


def avisar_uma_vez(msg: str):
    if msg not in _AVISOS:
        _AVISOS.add(msg)
        print(f"AVISO: {msg}", file=sys.stderr, flush=True)


# ----------------------------------------------------------------------------
# Geometria (tudo parte de WGS84; métrica só via geodésica ou plano local)
# ----------------------------------------------------------------------------
def area_geodesica_ha(geom) -> float:
    if geom.is_empty:
        return 0.0
    if hasattr(geom, "geoms"):
        return sum(area_geodesica_ha(g) for g in geom.geoms)
    return abs(GEOD.geometry_area_perimeter(geom)[0]) / 1e4


def crs_local(geom_ll) -> CRS:
    """Azimutal equidistante centrado no talhão: métrico, sem zona a escolher."""
    c = geom_ll.centroid
    return CRS.from_proj4(
        f"+proj=aeqd +lat_0={c.y:.7f} +lon_0={c.x:.7f} +x_0=0 +y_0=0 "
        f"+datum=WGS84 +units=m +no_defs"
    )


def func_transf(origem, destino):
    return Transformer.from_crs(origem, destino, always_xy=True).transform


def listar_arquivos(entradas) -> list[Path]:
    arquivos: list[Path] = []
    for e in entradas:
        p = Path(e)
        cands = sorted(p.rglob("*")) if p.is_dir() else [p]
        for c in cands:
            suf = c.suffix.lower()
            if suf == ".zip":
                d = tempfile.mkdtemp(prefix="zp_")
                extrair_zip_seguro(c, d)
                arquivos += [x for x in sorted(Path(d).rglob("*")) if x.suffix.lower() in EXTENSOES]
            elif suf in EXTENSOES:
                arquivos.append(c)
    return arquivos


def ler_campo(path: Path) -> gpd.GeoDataFrame:
    try:
        g = gpd.read_file(path)
    except Exception:
        g = gpd.read_file(path, driver="KML")
    if g.crs is None:
        avisar_uma_vez(f"{path.name}: sem CRS (.prj ausente); assumindo WGS84.")
        g = g.set_crs(4326)
    g = g.to_crs(4326)
    g = g[g.geometry.notna() & g.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    g["geometry"] = shapely.force_2d(g.geometry.values)
    inval = ~g.is_valid
    if inval.any():
        g.loc[inval, "geometry"] = g.loc[inval].geometry.make_valid()
    g = g.reset_index(drop=True)
    if "sample_id" in g.columns:
        g["sample_id"] = g["sample_id"].astype(str)
    else:
        g["sample_id"] = [f"c{i + 1}" for i in range(len(g))]
    return g[["sample_id", "geometry"]]


def preparar_geometrias(gdf, buffer_borda_m: float):
    """Devolve (crs_local, células_locais, união_local, células_internas_locais)."""
    ll = list(gdf.geometry.values)
    uniao_ll = unary_union(ll)
    crs_loc = crs_local(uniao_ll)
    f = func_transf("EPSG:4326", crs_loc)
    cel_loc = [shp_transform(f, g) for g in ll]
    uni_loc = unary_union(cel_loc)
    interior = uni_loc.buffer(-buffer_borda_m)
    if interior.is_empty:
        avisar_uma_vez("buffer de borda esvaziou o talhão; usando o talhão inteiro.")
        interior = uni_loc
    cel_int = [c.intersection(interior) for c in cel_loc]
    return crs_loc, cel_loc, uni_loc, cel_int


# ----------------------------------------------------------------------------
# Janelas de safra
# ----------------------------------------------------------------------------
def rotulo_janela(dt: datetime) -> str:
    if dt.month >= 10:
        return f"V{dt.year}"
    if dt.month <= 3:
        return f"V{dt.year - 1}"
    return f"I{dt.year}"


def limites_janela(rot: str):
    ano = int(rot[1:])
    if rot[0] == "V":
        return date(ano, 10, 1), date(ano + 1, 3, 31)
    return date(ano, 4, 1), date(ano, 9, 30)


def agrupar_por_janela(itens, ini: date, fim: date):
    grupos: dict[str, list] = {}
    for it in itens:
        if it.datetime is None:
            continue
        grupos.setdefault(rotulo_janela(it.datetime), []).append(it)
    out = []
    for rot in sorted(grupos, key=lambda r: limites_janela(r)[0]):
        a, b = limites_janela(rot)
        if a >= ini and b <= fim:
            out.append((rot, grupos[rot]))
    return out


# ----------------------------------------------------------------------------
# Sentinel-2 via STAC
# ----------------------------------------------------------------------------
def _cliente():
    from pystac_client import Client
    return Client.open(STAC_URL)


def buscar_itens(colecao, bbox, ini, fim, nuvem_max, tentativas=3):
    for t in range(tentativas):
        try:
            busca = _cliente().search(
                collections=[colecao], bbox=bbox,
                datetime=f"{ini}T00:00:00Z/{fim}T23:59:59Z",
                query={"eo:cloud_cover": {"lt": nuvem_max}},
            )
            return list(busca.items())
        except Exception:
            if t == tentativas - 1:
                raise
            time.sleep(5 * (t + 1))


def buscar_dem(bbox, colecao="cop-dem-glo-30", tentativas=3):
    for t in range(tentativas):
        try:
            return list(_cliente().search(collections=[colecao], bbox=bbox).items())
        except Exception:
            if t == tentativas - 1:
                raise
            time.sleep(5 * (t + 1))


def epsg_do_item(it):
    p = it.properties
    if p.get("proj:epsg"):
        return int(p["proj:epsg"])
    code = str(p.get("proj:code", ""))
    if code.upper().startswith("EPSG:"):
        return int(code.split(":")[1])
    return None


def epsg_dominante(itens):
    c = Counter(e for e in (epsg_do_item(i) for i in itens) if e)
    if not c:
        raise RuntimeError("Itens STAC sem informação de EPSG.")
    return c.most_common(1)[0][0]


def escala_offset(item, banda="red"):
    """(scale, offset) para converter DN em reflectância.
    Produtos com baseline >= 04.00 (processados desde 2022) têm offset -0.1."""
    a = item.assets.get(banda)
    rb = {}
    if a is not None:
        lst = a.extra_fields.get("raster:bands") or [{}]
        rb = lst[0] or {}
    if rb.get("scale") is not None and rb.get("offset") is not None:
        return float(rb["scale"]), float(rb["offset"])
    base = item.properties.get("s2:processing_baseline")
    try:
        return 1e-4, (-0.1 if float(base) >= 4.0 else 0.0)
    except (TypeError, ValueError):
        avisar_uma_vez("item sem raster:bands nem baseline; assumindo sem offset (conferir!).")
        return 1e-4, 0.0


def geobox_alinhado(bounds, epsg, res, margem):
    from odc.geo.geobox import GeoBox
    xmin, ymin, xmax, ymax = bounds
    bb = (math.floor((xmin - margem) / res) * res, math.floor((ymin - margem) / res) * res,
          math.ceil((xmax + margem) / res) * res, math.ceil((ymax + margem) / res) * res)
    return GeoBox.from_bbox(bb, crs=f"EPSG:{epsg}", resolution=res)


def carregar_ndvi(itens, gb):
    """NDVI (t, y, x) com NaN onde inválido. O offset de reflectância é tratado por item."""
    import odc.stac
    grupos: dict[tuple, list] = {}
    for it in itens:
        grupos.setdefault(escala_offset(it), []).append(it)
    fatias = []
    for (sc, off), its in grupos.items():
        ds = odc.stac.load(its, bands=["red", "nir", "scl"], geobox=gb, groupby="solar_day",
                           dtype="float32", resampling="nearest", fail_on_error=False)
        red = ds["red"].values * sc + off
        nir = ds["nir"].values * sc + off
        scl = ds["scl"].values
        den = nir + red
        with np.errstate(all="ignore"):
            ndvi = (nir - red) / den
        ok = np.isin(scl, SCL_VALIDO) & np.isfinite(ndvi) & (den > 0.02)
        fatias.append(np.where(ok, np.clip(ndvi, -1, 1), np.nan).astype("float32"))
    if not fatias:
        return np.empty((0,) + tuple(gb.shape), dtype="float32")
    return np.concatenate(fatias, axis=0)


def p90_janela(ndvi, obs_min):
    n = np.isfinite(ndvi).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = np.nanpercentile(ndvi, 90, axis=0) if ndvi.shape[0] else np.full(ndvi.shape[1:], np.nan)
    p = p.astype("float32")
    p[n < obs_min] = np.nan
    return p, n.astype("int16")


# ----------------------------------------------------------------------------
# Persistência
# ----------------------------------------------------------------------------
def metricas_persistencia(P, mascara, lim, ndvi_min, cob_min):
    """P: (janelas, y, x) com NDVI p90. Devolve métricas por pixel ou None."""
    total = max(int(mascara.sum()), 1)
    usadas, D = [], []
    for k in range(P.shape[0]):
        v = P[k]
        ok = mascara & np.isfinite(v)
        if ok.sum() / total < cob_min:
            continue
        mu = float(np.mean(v[ok]))
        if mu < ndvi_min:
            continue
        d = np.full(v.shape, np.nan, dtype="float32")
        d[ok] = (v[ok] - mu) / mu
        D.append(d)
        usadas.append(k)
    if not usadas:
        return None
    D = np.stack(D)
    n = np.isfinite(D).sum(axis=0)
    min_n = min(max(3, math.ceil(0.5 * len(usadas))), len(usadas))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mu_p, sd_p = np.nanmean(D, axis=0), np.nanstd(D, axis=0)
    den = np.maximum(n, 1)
    fneg = (np.sum(D < -lim, axis=0) / den).astype("float32")
    fpos = (np.sum(D > lim, axis=0) / den).astype("float32")
    inval = n < min_n
    for a in (mu_p, sd_p, fneg, fpos):
        a[inval] = np.nan
    return dict(usadas=usadas, D=D, mu=mu_p, sigma=sd_p, fneg=fneg, fpos=fpos, n=n)


# ----------------------------------------------------------------------------
# Relevo
# ----------------------------------------------------------------------------
_VIZ = [(-1, -1, math.sqrt(2)), (-1, 0, 1.0), (-1, 1, math.sqrt(2)), (0, -1, 1.0),
        (0, 1, 1.0), (1, -1, math.sqrt(2)), (1, 0, 1.0), (1, 1, math.sqrt(2))]


def preencher_depressoes(z, eps=1e-3):
    """Priority-flood com incremento eps (impõe gradiente em áreas planas)."""
    ny, nx = z.shape
    f = z.astype("float64").copy()
    vis = np.zeros(z.shape, bool)
    heap = []
    for i in range(ny):
        for j in (0, nx - 1):
            if not vis[i, j]:
                vis[i, j] = True
                heapq.heappush(heap, (f[i, j], i, j))
    for j in range(nx):
        for i in (0, ny - 1):
            if not vis[i, j]:
                vis[i, j] = True
                heapq.heappush(heap, (f[i, j], i, j))
    while heap:
        h, i, j = heapq.heappop(heap)
        for di, dj, _ in _VIZ:
            a, b = i + di, j + dj
            if 0 <= a < ny and 0 <= b < nx and not vis[a, b]:
                vis[a, b] = True
                f[a, b] = max(f[a, b], h + eps)
                heapq.heappush(heap, (f[a, b], a, b))
    return f


def acumulacao_d8(f):
    ny, nx = f.shape
    acc = np.ones((ny, nx))
    for flat in np.argsort(-f, axis=None):
        i, j = divmod(int(flat), nx)
        melhor, alvo = 0.0, None
        for di, dj, dist in _VIZ:
            a, b = i + di, j + dj
            if 0 <= a < ny and 0 <= b < nx:
                queda = (f[i, j] - f[a, b]) / dist
                if queda > melhor:
                    melhor, alvo = queda, (a, b)
        if alvo is not None:
            acc[alvo] += acc[i, j]
    return acc


def preencher_nan(z):
    if np.isfinite(z).all():
        return z
    if not np.isfinite(z).any():
        raise RuntimeError("DEM sem nenhum valor válido na área.")
    idx = ndi.distance_transform_edt(~np.isfinite(z), return_distances=False, return_indices=True)
    return z[tuple(idx)]


def derivar_relevo(z, res, raio_tpi_m=200, sigma_px=1.0):
    z = preencher_nan(np.asarray(z, dtype="float64"))
    zs = ndi.gaussian_filter(z, sigma_px, mode="nearest")
    gy, gx = np.gradient(zs, res)
    tanb = np.hypot(gx, gy)
    k = int(round(raio_tpi_m / res)) * 2 + 1
    tpi = zs - ndi.uniform_filter(zs, size=k, mode="nearest")
    acc = acumulacao_d8(preencher_depressoes(zs))
    a = acc * res                                  # área específica de contribuição (m)
    twi = np.log(a / np.maximum(tanb, 1e-3))
    return dict(z=zs, decliv_pct=100 * tanb, tpi200=tpi, twi=twi)


def carregar_dem(bbox, gb):
    import odc.stac
    itens = buscar_dem(bbox)
    if not itens:
        raise RuntimeError("Nenhum tile do Copernicus DEM encontrado para a área.")
    chave = "data" if "data" in itens[0].assets else next(
        k for k, a in itens[0].assets.items() if "data" in (a.roles or []))
    ds = odc.stac.load(itens, bands=[chave], geobox=gb, groupby="id", dtype="float32",
                       resampling="bilinear", fail_on_error=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmax(ds[chave].values, axis=0)


# ----------------------------------------------------------------------------
# Agregação por célula
# ----------------------------------------------------------------------------
def rotular_celulas(polys, affine, shape):
    formas = [(g, i + 1) for i, g in enumerate(polys) if not g.is_empty]
    if not formas:
        return np.zeros(shape, dtype="int32")
    return rasterio.features.rasterize(formas, out_shape=shape, transform=affine,
                                       fill=0, dtype="int32", all_touched=False)


def agregar(rot, feats, n_cel, sd_de=(), sufixo_npix="n_pix"):
    linhas = []
    for i in range(1, n_cel + 1):
        m = rot == i
        r = {sufixo_npix: int(m.sum())}
        for nome, arr in feats.items():
            v = arr[m]
            v = v[np.isfinite(v)]
            r[nome] = float(v.mean()) if v.size else np.nan
            if nome in sd_de:
                r[nome + "_sd_intra"] = float(v.std()) if v.size > 1 else np.nan
        linhas.append(r)
    return linhas


# ----------------------------------------------------------------------------
# Pipeline por talhão
# ----------------------------------------------------------------------------
def slug(s):
    return "".join(c if c.isalnum() else "_" for c in s)[:60]


def processar_campo(campo, gdf, args):
    t0 = time.time()
    ini, fim = date.fromisoformat(args.inicio), date.fromisoformat(args.fim)
    crs_loc, cel_loc, uni_loc, cel_int = preparar_geometrias(gdf, args.buffer_borda)
    uniao_ll = unary_union(list(gdf.geometry.values))
    bbox = tuple(uniao_ll.buffer(0.001).bounds)

    chave = hashlib.md5(f"{campo}|{bbox}|{ini}|{fim}|{args.colecao}|{args.nuvem_max}|"
                        f"{args.obs_min}|{len(gdf)}".encode()).hexdigest()[:10]
    cache = Path(args.cache) / f"{slug(campo)}_{chave}.npz"

    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        P, labels, aff6, epsg = z["P"], list(z["labels"]), tuple(z["aff"]), int(z["epsg"])
        log(f"  cache: {cache.name} ({len(labels)} janelas)")
        from affine import Affine
        aff = Affine(*aff6)
    else:
        itens = buscar_itens(args.colecao, bbox, ini, fim, args.nuvem_max)
        if not itens:
            raise RuntimeError("Nenhuma cena Sentinel-2 encontrada.")
        epsg = epsg_dominante(itens)
        itens = [i for i in itens if epsg_do_item(i) == epsg]
        log(f"  {len(itens)} cenas (EPSG:{epsg}), primeira em "
            f"{min(i.datetime for i in itens).date()}")
        f_r = func_transf(crs_loc, f"EPSG:{epsg}")
        uni_r = shp_transform(f_r, uni_loc)
        gb = geobox_alinhado(uni_r.bounds, epsg, RES_S2, margem=40)
        aff = gb.affine
        Ps, labels = [], []
        for rot, its in agrupar_por_janela(itens, ini, fim):
            ndvi = carregar_ndvi(its, gb)
            p, n = p90_janela(ndvi, args.obs_min)
            Ps.append(p)
            labels.append(rot)
            log(f"    {rot}: {len(its)} cenas, obs válidas/pixel (mediana) = "
                f"{int(np.median(n))}")
        if not Ps:
            raise RuntimeError("Nenhuma janela completa dentro do período.")
        P = np.stack(Ps)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, P=P, labels=np.array(labels),
                            aff=np.array(tuple(aff)[:6]), epsg=epsg)

    f_r = func_transf(crs_loc, f"EPSG:{epsg}")
    cel_int_r = [shp_transform(f_r, c) for c in cel_int]
    uni_r = shp_transform(f_r, uni_loc)
    ny, nx = P.shape[1:]
    rot = rotular_celulas(cel_int_r, aff, (ny, nx))
    mascara = rot > 0

    met = metricas_persistencia(P, mascara, args.limiar_anomalia, args.ndvi_min_janela,
                                args.cobertura_min)
    if met is None:
        raise RuntimeError("Nenhuma janela com lavoura/cobertura suficiente.")
    log(f"  janelas usadas: {len(met['usadas'])}/{len(labels)}")
    feats_s2 = {"ndvi_mu": met["mu"], "ndvi_sigma": met["sigma"],
                "ndvi_fneg": met["fneg"], "ndvi_fpos": met["fpos"]}
    n_cel = len(gdf)
    s2 = agregar(rot, feats_s2, n_cel, sd_de=("ndvi_mu",))

    # relevo
    bbox_dem = tuple(uniao_ll.buffer(0.02).bounds)       # ~2 km de margem para a acumulação
    gb_dem = geobox_alinhado(uni_r.bounds, epsg, RES_DEM, margem=1500)
    z = carregar_dem(bbox_dem, gb_dem)
    rel = derivar_relevo(z, RES_DEM)
    rot_d = rotular_celulas(cel_int_r, gb_dem.affine, z.shape)
    rel["elev_rel"] = rel["z"] - float(np.mean(rel["z"][rot_d > 0]))
    dem = agregar(rot_d, {k: rel[k] for k in ("decliv_pct", "tpi200", "elev_rel", "twi")},
                  n_cel, sd_de=("elev_rel", "twi"), sufixo_npix="n_pix_dem")

    # tabelas
    linhas, jan = [], []
    for i in range(n_cel):
        c = gdf.geometry.iloc[i]
        r = {"campo": campo, "celula": i + 1, "sample_id": gdf["sample_id"].iloc[i],
             "lon": c.centroid.x, "lat": c.centroid.y, "area_ha": area_geodesica_ha(c),
             "n_janelas": len(met["usadas"])}
        r.update(s2[i])
        r.update(dem[i])
        linhas.append(r)
        m = rot == i + 1
        for j, k in enumerate(met["usadas"]):
            jan.append({"campo": campo, "celula": i + 1, "sample_id": gdf["sample_id"].iloc[i],
                        "janela": labels[k],
                        "ndvi_p90": float(np.nanmean(P[k][m])) if m.any() else np.nan,
                        "anomalia_rel": float(np.nanmean(met["D"][j][m])) if m.any() else np.nan})
    log(f"  ok em {time.time() - t0:.0f}s")
    return pd.DataFrame(linhas), pd.DataFrame(jan)


# ----------------------------------------------------------------------------
# Autoteste (sem rede)
# ----------------------------------------------------------------------------
def _selftest():
    from shapely.geometry import Polygon, box
    ok = True

    def check(nome, cond, extra=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"[{'OK' if cond else 'FALHA'}] {nome} {extra}")

    # 1. geometria: quadrado 1 km dividido em 4 células, em Mato Grosso
    centro = (-54.85, -10.39)
    crs_l = crs_local(shapely.Point(*centro))
    f_back = func_transf(crs_l, "EPSG:4326")
    cel = [shp_transform(f_back, box(x, y, x + 500, y + 500))
           for x in (-500, 0) for y in (-500, 0)]
    gdf = gpd.GeoDataFrame({"sample_id": list("abcd")}, geometry=cel, crs=4326)
    a = [area_geodesica_ha(g) for g in gdf.geometry]
    check("área geodésica das células (esperado ~25 ha)", all(abs(x - 25) < 0.2 for x in a),
          f"{[round(x, 2) for x in a]}")
    _, _, uni, cint = preparar_geometrias(gdf, 20)
    area_int = sum(c.area for c in cint) / 1e4
    check("buffer de borda 20 m (esperado ~ 96^2/100 = 92,2 ha)", abs(area_int - 92.16) < 0.5,
          f"{area_int:.2f} ha")

    # 2. persistência: mancha 20% abaixo, persistente em todas as janelas.
    #    A mancha (25% da área) também puxa a média do talhão: mu esperado = -0,158.
    rng = np.random.default_rng(1)
    ny = nx = 100
    P = []
    for _ in range(6):
        v = 0.80 + rng.normal(0, 0.015, (ny, nx))
        v[:50, :50] *= 0.80
        P.append(v.astype("float32"))
    P.append(np.full((ny, nx), 0.20, dtype="float32"))           # pousio: deve ser descartada
    P = np.stack(P)
    m = np.ones((ny, nx), bool)
    r = metricas_persistencia(P, m, 0.10, 0.5, 0.7)
    check("janela de pousio descartada", len(r["usadas"]) == 6, f"usadas={len(r['usadas'])}")
    mu_m, fn_m = np.nanmean(r["mu"][:50, :50]), np.nanmean(r["fneg"][:50, :50])
    check("mancha persistente: mu ~ -0,158 e fneg ~ 1", abs(mu_m + 0.158) < 0.01 and fn_m > 0.95,
          f"mu={mu_m:.3f} fneg={fn_m:.2f}")
    check("fora da mancha: fneg ~ 0", np.nanmean(r["fneg"][60:, 60:]) < 0.05)

    # 3. relevo: vale em V com declive suave; TWI no eixo do vale > nas encostas
    y, x = np.mgrid[0:80, 0:80].astype(float)
    z = 100 + 0.8 * np.abs(x - 40) + 0.2 * y
    rel = derivar_relevo(z, 30)
    check("TWI maior no eixo do vale que na encosta",
          rel["twi"][40, 40] > rel["twi"][40, 10] + 1.0,
          f"vale={rel['twi'][40, 40]:.2f} encosta={rel['twi'][40, 10]:.2f}")
    acc = acumulacao_d8(preencher_depressoes(ndi.gaussian_filter(z, 1, mode="nearest")))
    check("acumulação conserva área (máx ~ nº de células)", acc.max() <= z.size + 1e-6)

    # 4. offset de reflectância
    import pystac
    def fake(base, rb=None):
        it = pystac.Item("x", None, None, datetime(2022, 3, 1), {"s2:processing_baseline": base})
        extra = {"raster:bands": [rb]} if rb else {}
        it.add_asset("red", pystac.Asset("x.tif", extra_fields=extra))
        return it
    check("offset por baseline 05.11 = -0.1", escala_offset(fake("05.11")) == (1e-4, -0.1))
    check("offset por baseline 03.00 = 0", escala_offset(fake("03.00")) == (1e-4, 0.0))
    check("offset lido de raster:bands", escala_offset(fake("05.11", {"scale": 1e-4, "offset": -0.1})) == (1e-4, -0.1))

    # 5. leitura de imagens com itens STAC locais (offset + SCL + nodata)
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import mapping
    epsg, x0, y0, n = 31981, 600000.0, 8850000.0, 200
    tmp = Path(tempfile.mkdtemp(prefix="zp_test_"))
    t_ll = Transformer.from_crs(epsg, 4326, always_xy=True)
    poly = Polygon([t_ll.transform(a, b) for a, b in
                    [(x0, y0), (x0 + n * 10, y0), (x0 + n * 10, y0 - n * 10), (x0, y0 - n * 10)]])

    def tif(path, arr, dtype, transform, crs, nodata=0):
        with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                           dtype=dtype, crs=crs, transform=transform, nodata=nodata) as d:
            d.write(arr.astype(dtype), 1)

    itens = []
    for k, (dtm, base) in enumerate([(datetime(2021, 12, 1, 13), "03.00"),
                                     (datetime(2022, 3, 1, 13), "05.11")]):
        off = -0.1 if float(base) >= 4 else 0.0
        red = np.full((n, n), (0.05 - off) / 1e-4)
        nir = np.full((n, n), (0.40 - off) / 1e-4)
        scl = np.full((n, n), 4)
        scl[10:50, 10:50] = 9                                     # nuvem
        for nm, arr, dt_ in [("red", red, "uint16"), ("nir", nir, "uint16"), ("scl", scl, "uint8")]:
            tif(tmp / f"{k}_{nm}.tif", arr, dt_, from_origin(x0, y0, 10, 10), f"EPSG:{epsg}")
        it = pystac.Item(f"s{k}", mapping(poly), list(poly.bounds), dtm,
                         {"s2:processing_baseline": base, "proj:epsg": epsg})
        for nm, dt_ in [("red", "uint16"), ("nir", "uint16"), ("scl", "uint8")]:
            rb = {"nodata": 0, "data_type": dt_}
            if nm != "scl":
                rb.update(scale=1e-4, offset=off)
            it.add_asset(nm, pystac.Asset(str(tmp / f"{k}_{nm}.tif"), roles=["data"],
                         media_type="image/tiff; application=geotiff",
                         extra_fields={"proj:shape": [n, n], "proj:epsg": epsg,
                                       "proj:transform": [10, 0, x0, 0, -10, y0, 0, 0, 1],
                                       "raster:bands": [rb]}))
        itens.append(it)
    gb = geobox_alinhado((x0 + 100, y0 - 1500, x0 + 1500, y0 - 100), epsg, 10, 0)
    nd = carregar_ndvi(itens, gb)
    v = nd[np.isfinite(nd)]
    check("NDVI idêntico nas duas baselines (0,778)", v.size > 0 and np.allclose(v, (0.40 - 0.05) / (0.40 + 0.05), atol=1e-3),
          f"min={v.min():.3f} max={v.max():.3f}")
    check("pixels com nuvem (SCL 9) mascarados", np.isnan(nd[:, 5:35, 5:35]).all())

    # 6. DEM sintético em 4326 reamostrado para a grade UTM
    import odc.stac
    lon0, lat0, d = poly.bounds[0] - 0.01, poly.bounds[3] + 0.01, 1 / 3600
    nn = 1200
    yy, xx = np.mgrid[0:nn, 0:nn]
    dem = (300 + 0.05 * xx).astype("float32")
    tif(tmp / "dem.tif", dem, "float32", from_origin(lon0, lat0, d, d), "EPSG:4326", nodata=-9999)
    bb = (lon0, lat0 - nn * d, lon0 + nn * d, lat0)
    dit = pystac.Item("dem", mapping(box(*bb)), list(bb), datetime(2021, 4, 22),
                      {"proj:epsg": 4326})
    dit.add_asset("data", pystac.Asset(str(tmp / "dem.tif"), roles=["data"],
                  media_type="image/tiff; application=geotiff",
                  extra_fields={"proj:shape": [nn, nn], "proj:epsg": 4326,
                                "proj:transform": [d, 0, lon0, 0, -d, lat0, 0, 0, 1],
                                "raster:bands": [{"nodata": -9999, "data_type": "float32"}]}))
    ds = odc.stac.load([dit], bands=["data"], geobox=geobox_alinhado(
        (x0 + 100, y0 - 1500, x0 + 1500, y0 - 100), epsg, 30, 0), groupby="id",
        dtype="float32", resampling="bilinear")
    zz = np.nanmax(ds["data"].values, axis=0)
    check("DEM reamostrado para UTM (valores plausíveis, sem NaN)",
          np.isfinite(zz).all() and 300 < zz.min() < zz.max() < 400, f"{zz.min():.0f}-{zz.max():.0f} m")

    print("\nTODOS OS TESTES PASSARAM" if ok else "\nHÁ TESTES COM FALHA")
    return 0 if ok else 1


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entradas", nargs="*", help="pastas, .zip, .shp, .kml, .geojson (um arquivo = um talhão)")
    ap.add_argument("--saida", default="celulas_features.csv")
    ap.add_argument("--inicio", default="2017-01-01")
    ap.add_argument("--fim", default=date.today().isoformat())
    ap.add_argument("--buffer-borda", type=float, default=20.0, help="m, recuo da borda do talhão")
    ap.add_argument("--colecao", default="sentinel-2-l2a")
    ap.add_argument("--nuvem-max", type=float, default=80.0, help="%% de nuvem máx. por cena")
    ap.add_argument("--obs-min", type=int, default=3, help="obs. válidas mínimas por pixel e janela")
    ap.add_argument("--limiar-anomalia", type=float, default=0.10)
    ap.add_argument("--ndvi-min-janela", type=float, default=0.5)
    ap.add_argument("--cobertura-min", type=float, default=0.7)
    ap.add_argument("--cache", default="cache_zp")
    ap.add_argument("--so-campo", default=None, help="processa só talhões cujo nome contém este texto")
    ap.add_argument("--so-geometria", action="store_true", help="valida entradas, sem rede")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(_selftest())
    if not args.entradas:
        ap.error("informe ao menos uma entrada (ou use --selftest)")

    arquivos = listar_arquivos(args.entradas)
    if args.so_campo:
        arquivos = [a for a in arquivos if args.so_campo.lower() in a.stem.lower()]
    if not arquivos:
        sys.exit("Nenhum arquivo de geometria encontrado.")

    todos, janelas, erros = [], [], []
    for path in arquivos:
        campo = path.stem
        log(f"\n== {campo}")
        try:
            gdf = ler_campo(path)
            if args.so_geometria:
                _, _, uni, cint = preparar_geometrias(gdf, args.buffer_borda)
                a_tot = sum(area_geodesica_ha(g) for g in gdf.geometry)
                a_int = sum(c.area for c in cint) / 1e4
                log(f"  {len(gdf)} células | {a_tot:.1f} ha | após recuo de {args.buffer_borda:g} m: {a_int:.1f} ha")
                continue
            df, jn = processar_campo(campo, gdf, args)
            todos.append(df)
            janelas.append(jn)
        except Exception as e:
            erros.append((campo, repr(e)))
            log(f"  ERRO: {e!r}")

    if args.so_geometria:
        return
    if todos:
        saida = Path(args.saida)
        pd.concat(todos, ignore_index=True).to_csv(saida, index=False)
        pd.concat(janelas, ignore_index=True).to_csv(saida.with_name(saida.stem + "_janelas.csv"), index=False)
        log(f"\nGravado: {saida} e {saida.stem}_janelas.csv")
    for campo, e in erros:
        log(f"FALHOU: {campo}: {e}")
    if erros:
        sys.exit(1)


if __name__ == "__main__":
    main()
