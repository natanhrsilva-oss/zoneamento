#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zonas_core.py — motor do zoneamento direcionado para amostragem TruBio
=====================================================================
Entrada : limite do talhão (KML/KMZ, shapefile, GeoJSON), em WGS84.
Saída   : KML com as zonas propostas (uma amostra composta por zona), mais camadas de
          estabilidade e de interesse (queda de vigor), RGB, painel e tabela.

Fontes públicas e gratuitas (Earth Search/AWS via STAC): Sentinel-2 L2A e Copernicus DEM.

COMO O TRÁFEGO DE DADOS É MINIMIZADO
  1. Busca STAC só com metadados (campos pesados excluídos) e UMA busca por fazenda.
  2. Só são lidas as bandas necessárias: vermelho + infravermelho próximo (NDVI) e a
     banda SCL. Não se baixa cena inteira: cada leitura é um recorte do talhão
     (requisições HTTP por faixa de bytes, formato COG).
  3. Triagem de nuvem ANTES do NDVI: lê-se só a SCL (20 m, a banda mais leve) do recorte.
     Cenas com o talhão encoberto nunca têm vermelho/NIR baixados.
  4. Orçamento por janela: para quando já há `obs_por_janela` cenas boas (padrão 6),
     começando pelas de menor nuvem. Não se baixa a série inteira.
  5. Resolução de trabalho de 20 m (padrão): lê o nível de overview do COG.
  6. Cache em disco por cena e por talhão: reanálises e novas safras só baixam o que
     ainda não existe. A RGB é UMA cena (a mais limpa e recente).
  7. Sem mapas-base externos no app (tiles=None): a imagem de fundo é a própria RGB.

CONVENÇÕES: tudo em WGS84 para o usuário; áreas geodésicas; operações métricas no CRS
nativo da imagem (escolhido a partir da própria cena, nunca pelo usuário).
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import io
import json
import math
import re
import sys
import tempfile
import threading
import time
import warnings
import zipfile
from dataclasses import dataclass, field, replace
from datetime import date, datetime, timedelta
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd

try:
    import geopandas as gpd
    import rasterio
    import rasterio.features
    import rasterio.warp
    import shapely
    from pyproj import Transformer
    from rasterio.enums import Resampling
    from rasterio.transform import Affine
    from rasterio.windows import Window
    from scipy import ndimage as ndi
    from shapely.geometry import MultiPolygon, Polygon, mapping, shape
    from shapely.geometry.polygon import orient
    from shapely.ops import transform as shp_transform
    from shapely.ops import unary_union
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.feature_extraction.image import grid_to_graph
except ImportError as e:  # pragma: no cover
    sys.exit(f"Dependência ausente: {e}. Veja requirements_zonas.txt")

import zonas_persistencia as zp
from arquivos_seguros import extrair_zip_seguro

RES_DEM = 30
_RIO_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF",
    GDAL_HTTP_MULTIPLEX="YES", GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", VSI_CACHE="TRUE",
    VSI_CACHE_SIZE="25000000", GDAL_HTTP_MAX_RETRY="3", GDAL_HTTP_RETRY_DELAY="2",
    AWS_NO_SIGN_REQUEST="YES",
)
_EXCLUIR_ASSETS = [f"assets.{a}" for a in (
    "coastal", "rededge1", "rededge2", "rededge3", "nir08", "nir09", "swir16", "swir22", "aot", "wvp",
    "thumbnail", "granule_metadata", "product_metadata", "tileinfo_metadata", "preview",
    "visual", "red-jp2", "green-jp2", "blue-jp2", "nir-jp2", "scl-jp2", "coastal-jp2",
    "rededge1-jp2", "rededge2-jp2", "rededge3-jp2", "nir08-jp2", "nir09-jp2", "swir16-jp2",
    "swir22-jp2", "aot-jp2", "wvp-jp2", "visual-jp2")]


# =============================================================================
# Parâmetros, estruturas
# =============================================================================
@dataclass
class Params:
    densidade_ha: float = 12.0        # 1 amostra a cada X ha (alvo)
    area_min_ha: float = 6.0          # limites FLEXÍVEIS (penalizados, não rígidos)
    area_max_ha: float = 20.0
    rigidez_area: float = 0.05        # 0 = ignora os limites; 0,05 = padrão; 0,2 = bem rígido
    piso_ha: float = 2.0              # abaixo disso a zona é fundida (inviável de amostrar)
    inicio: str = "2017-01-01"
    fim: str = field(default_factory=lambda: date.today().isoformat())
    res_m: int = 20                   # 20 (padrão, menos tráfego) ou 10
    obs_por_janela: int = 6           # orçamento de cenas boas por janela de safra
    obs_min_janela: int = 3
    obs_min_pixel: int = 3
    nuvem_max_cena: float = 60.0      # filtro de metadados (sem baixar nada)
    claro_min: float = 0.6            # fração do talhão sem nuvem para aceitar a cena
    buffer_borda_m: float = 20.0
    limiar_anomalia: float = 0.10     # desvio relativo à mediana do talhão
    ndvi_min_janela: float = 0.5      # janelas sem lavoura são descartadas
    cobertura_min: float = 0.7
    base: str = "mediana"             # "mediana" ou "media"
    f_interesse: float = 0.5          # queda persistente: abaixo do limiar em >= 50% das janelas
    f_estavel: float = 0.15           # estável: desvio (+/-) em <= 15% das janelas
    n_recente: int = 3                # queda recente: abaixo do limiar em >= 2/3 das últimas n janelas e normal antes
    peso_relevo: float = 0.5
    compacidade: float = 0.5          # 0 = zonas seguem só as covariáveis; 1 = formas bem compactas
    suavizacao: float = 0.8           # 0 = contornos em degraus de pixel; 1 = contornos bem arredondados
    area_min_camada_ha: float = 1.0
    dias_rgb: int = 120
    cache_dir: str = "cache_zonas"
    workers: int = 6


@dataclass
class Cena:
    id: str
    datetime: datetime
    nuvem: float
    epsg: int
    item: object = None


@dataclass(frozen=True)
class Grade:
    epsg: int
    x0: float
    y0: float
    res: float
    ny: int
    nx: int

    @property
    def affine(self):
        return Affine(self.res, 0, self.x0, 0, -self.res, self.y0)

    @property
    def shape(self):
        return (self.ny, self.nx)

    def bounds(self):
        return (self.x0, self.y0 - self.ny * self.res, self.x0 + self.nx * self.res, self.y0)

    def geobox(self):
        from odc.geo.geobox import GeoBox
        return GeoBox.from_bbox(self.bounds(), crs=f"EPSG:{self.epsg}", resolution=self.res)


def grade_alinhada(bounds, epsg, res, margem, alinhar=None):
    a = alinhar or res
    xmin, ymin, xmax, ymax = bounds
    x0 = math.floor((xmin - margem) / a) * a
    x1 = math.ceil((xmax + margem) / a) * a
    y1 = math.floor((ymin - margem) / a) * a
    y0 = math.ceil((ymax + margem) / a) * a
    return Grade(epsg, float(x0), float(y0), float(res), int(round((y0 - y1) / res)),
                 int(round((x1 - x0) / res)))


@dataclass
class Trafego:
    cenas_encontradas: int = 0
    triagens_rede: int = 0
    triagens_cache: int = 0
    ndvi_rede: int = 0
    ndvi_cache: int = 0
    cenas_descartadas_nuvem: int = 0
    rgb_rede: int = 0
    rgb_cache: int = 0
    dem_rede: int = 0
    bytes_medidos: float | None = None

    @property
    def leituras_banda(self):
        """Leituras de recorte feitas NESTA execução pela rede (cena x banda): SCL=1, NDVI=2, RGB=3."""
        return self.triagens_rede + 2 * self.ndvi_rede + 3 * self.rgb_rede

    @property
    def leituras_necessarias(self):
        """Leituras que a configuração exige, contando as que já estavam no cache."""
        return (self.triagens_rede + self.triagens_cache) + 2 * (self.ndvi_rede + self.ndvi_cache) + 3 * (self.rgb_rede + self.rgb_cache)

    def somar(self, o):
        for k in ("cenas_encontradas", "triagens_rede", "triagens_cache", "ndvi_rede", "ndvi_cache",
                  "cenas_descartadas_nuvem", "rgb_rede", "rgb_cache", "dem_rede"):
            setattr(self, k, getattr(self, k) + getattr(o, k))
        if o.bytes_medidos is not None:
            self.bytes_medidos = (self.bytes_medidos or 0) + o.bytes_medidos


class MedidorRede:
    """Mede bytes recebidos pela máquina durante o bloco (aproximado: inclui outros processos)."""
    def __enter__(self):
        try:
            import psutil
            self._ps = psutil
            self._b0 = psutil.net_io_counters().bytes_recv
        except Exception:
            self._ps = None
        self.bytes = None
        return self

    def __exit__(self, *a):
        if self._ps is not None:
            self.bytes = float(self._ps.net_io_counters().bytes_recv - self._b0)


# =============================================================================
# Leitura de recortes COG (só o necessário)
# =============================================================================
def ler_janela(href, grade, resamp=Resampling.average, overview_level=None):
    """Lê apenas o recorte da grade. Com overview_level abre o nível reduzido do COG
    (menos blocos baixados); se o arquivo não tiver overview, cai para decimação."""
    with rasterio.Env(**_RIO_ENV):
        src = None
        if overview_level is not None:
            try:
                src = rasterio.open(href, overview_level=overview_level)
            except rasterio.errors.RasterioIOError:
                src = None
        if src is None:
            src = rasterio.open(href)
        with src:
            l, b, r, t = grade.bounds()
            col = round((l - src.transform.c) / src.transform.a)
            row = round((t - src.transform.f) / src.transform.e)
            w = max(1, round((r - l) / src.transform.a))
            h = max(1, round((t - b) / -src.transform.e))
            return src.read(1, window=Window(col, row, w, h), out_shape=(grade.ny, grade.nx),
                            resampling=resamp, boundless=True, fill_value=0)


def ndvi_de(red_dn, nir_dn, scl, sc, off):
    red = red_dn.astype("float32") * sc + off
    nir = nir_dn.astype("float32") * sc + off
    den = nir + red
    with np.errstate(all="ignore"):
        v = (nir - red) / den
    ok = np.isin(scl, zp.SCL_VALIDO) & (red_dn > 0) & (nir_dn > 0) & (den > 0.02)
    return np.where(ok, np.clip(v, -1, 1), np.nan).astype("float32")


# =============================================================================
# Fontes de dados
# =============================================================================
class FonteSTAC:
    nome = "Earth Search (AWS) · STAC"
    sintetica = False

    def __init__(self, colecao="sentinel-2-l2a", url=zp.STAC_URL):
        self.colecao, self.url = colecao, url
        self._usar_fields = True

    def buscar(self, bbox, ini, fim, nuvem_max):
        from pystac_client import Client
        ultimo = None
        for tent in range(4):
            try:
                cli = Client.open(self.url)
                kw = dict(collections=[self.colecao], bbox=bbox,
                          datetime=f"{ini}T00:00:00Z/{fim}T23:59:59Z",
                          query={"eo:cloud_cover": {"lt": nuvem_max}})
                if self._usar_fields:
                    kw["fields"] = {"exclude": _EXCLUIR_ASSETS}
                itens = list(cli.search(**kw).items())
                break
            except Exception as e:
                ultimo = e
                self._usar_fields = False      # servidor pode não aceitar 'fields'
                time.sleep(3 * (tent + 1))
        else:
            raise RuntimeError(f"Falha na busca STAC: {ultimo!r}")
        cenas = []
        for it in itens:
            e = zp.epsg_do_item(it)
            if it.datetime is None or not e:
                continue
            cenas.append(Cena(it.id, it.datetime, float(it.properties.get("eo:cloud_cover", 100.0)), e, it))
        return cenas

    def _href(self, cena, banda):
        return cena.item.assets[banda].href

    def claro(self, cena, grade20, mask20):
        scl = ler_janela(self._href(cena, "scl"), grade20, Resampling.nearest)
        return float(np.isin(scl, zp.SCL_VALIDO)[mask20].mean()) if mask20.any() else 0.0

    def ndvi(self, cena, grade):
        sc, off = zp.escala_offset(cena.item)
        ov = 0 if grade.res >= 20 else None
        red = ler_janela(self._href(cena, "red"), grade, Resampling.average, ov)
        nir = ler_janela(self._href(cena, "nir"), grade, Resampling.average, ov)
        scl = ler_janela(self._href(cena, "scl"), grade, Resampling.nearest)
        return ndvi_de(red, nir, scl, sc, off)

    def rgb(self, cena, grade):
        sc, off = zp.escala_offset(cena.item)
        bandas = [ler_janela(self._href(cena, b), grade, Resampling.average) for b in ("red", "green", "blue")]
        refl = np.stack([b.astype("float32") * sc + off for b in bandas])
        validos = np.all(np.stack(bandas) > 0, axis=0)
        img = np.clip(refl / 0.28, 0, 1) ** (1 / 1.8)
        img[:, ~validos] = 0
        return (img * 255).astype("uint8")

    def dem(self, bbox, grade):
        return zp.carregar_dem(bbox, grade.geobox())


class FonteSintetica:
    """Dados sintéticos para demonstração e testes (sem rede). Padrões em coordenadas
    métricas, portanto idênticos em qualquer resolução."""
    nome = "SINTÉTICA (demonstração)"
    sintetica = True

    def __init__(self, semente=11):
        self.semente = semente

    @staticmethod
    def _seed(txt):
        return int(hashlib.md5(txt.encode()).hexdigest()[:8], 16)

    def buscar(self, bbox, ini, fim, nuvem_max):
        lon, lat = (bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2
        epsg = 31960 + int((lon + 180) // 6) + 1
        rng = np.random.default_rng(self.semente)
        d, fim_d, cenas = date.fromisoformat(ini), date.fromisoformat(fim), []
        while d <= fim_d:
            seco = d.month in (5, 6, 7, 8, 9)
            nuv = float(np.clip(rng.gamma(2.0, 4.0 if seco else 22.0), 0, 100))
            if nuv < nuvem_max:
                cenas.append(Cena(f"SYN_{d:%Y%m%d}", datetime(d.year, d.month, d.day, 13, 30), nuv, epsg))
            d += timedelta(days=5)
        return cenas

    def _eixos(self, grade):
        X = grade.x0 + grade.res * (np.arange(grade.nx) + 0.5)
        Y = grade.y0 - grade.res * (np.arange(grade.ny) + 0.5)
        X, Y = np.meshgrid(X, Y)
        l, b, r, t = grade.bounds()
        return X, Y, (l + r) / 2, (b + t) / 2, (r - l) / 2, (t - b) / 2

    def _padrao(self, grade):
        X, Y, cx, cy, hw, hh = self._eixos(grade)
        s = 0.5 * np.sin((X - cx) / 330 + 0.7) * np.cos((Y - cy) / 410 - 0.3) + 0.35 * np.sin(((X - cx) + (Y - cy)) / 260)
        rr = 0.28 * min(hw, hh)
        pers = np.exp(-(((X - cx + 0.35 * hw) ** 2 + (Y - cy - 0.2 * hh) ** 2) / (2 * rr ** 2)))
        rec = np.exp(-(((X - cx - 0.3 * hw) ** 2 + (Y - cy + 0.3 * hh) ** 2) / (2 * (0.8 * rr) ** 2)))
        return s, pers, rec

    def _nuvem(self, cena, grade):
        rng = np.random.default_rng(self._seed(cena.id))
        frac = np.clip(cena.nuvem / 100 * rng.uniform(0.7, 1.3), 0, 1)
        r = ndi.gaussian_filter(rng.normal(size=grade.shape), 6 / (grade.res / 10))
        return r > np.quantile(r, 1 - frac) if frac > 0.01 else np.zeros(grade.shape, bool)

    def claro(self, cena, grade20, mask20):
        return float((~self._nuvem(cena, grade20))[mask20].mean())

    def _ndvi_limpo(self, cena, grade):
        s, pers, rec = self._padrao(grade)
        base = 0.32 if (zp.rotulo_janela(cena.datetime)[0] == "I" and int(zp.rotulo_janela(cena.datetime)[1:]) % 2 == 0) else 0.80
        recente = 1.0 if cena.datetime >= datetime(2025, 4, 1) else 0.0
        return base + (0.045 * s - 0.20 * pers - 0.28 * rec * recente) * (1 if base > 0.5 else 0.2)

    def ndvi(self, cena, grade):
        rng = np.random.default_rng(self._seed(cena.id + "n"))
        v = self._ndvi_limpo(cena, grade) + rng.normal(0, 0.02, grade.shape)
        v[self._nuvem(cena, grade)] = np.nan
        return v.astype("float32")

    def rgb(self, cena, grade):
        v = self._ndvi_limpo(replace(cena, datetime=datetime(2026, 1, 15)), grade)
        f = np.clip((v - 0.3) / 0.5, 0, 1)
        rng = np.random.default_rng(self._seed(cena.id + "r"))
        tex = rng.normal(0, 4, grade.shape)
        img = np.stack([150 - 110 * f + tex, 115 + 20 * f + tex, 85 - 45 * f + tex])
        return np.clip(img, 0, 255).astype("uint8")

    def dem(self, bbox, grade):
        X, Y, cx, cy, hw, hh = self._eixos(grade)
        return (380 + 0.004 * (X - cx) + 0.002 * (Y - cy) + 7 * np.exp(-(((X - cx - 300) ** 2 + (Y - cy) ** 2) / (2 * 400 ** 2)))
                - 5 * np.exp(-(((X - cx + 200) ** 2 + (Y - cy - 250) ** 2) / (2 * 300 ** 2)))).astype("float32")


# =============================================================================
# Cache em disco
# =============================================================================
def _slug(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s))[:80]


class Cache:
    def __init__(self, raiz, chave):
        self.dir = Path(raiz) / chave
        self.dir.mkdir(parents=True, exist_ok=True)
        self._tri_path = self.dir / "triagem.json"
        self.tri = json.loads(self._tri_path.read_text()) if self._tri_path.exists() else {}
        self._lock = threading.Lock()

    def set_tri(self, cid, v):
        with self._lock:
            self.tri[cid] = float(v)

    def salvar_tri(self):
        with self._lock:
            self._tri_path.write_text(json.dumps(self.tri))

    def get_ndvi(self, cid):
        p = self.dir / f"n_{_slug(cid)}.npy"
        return np.load(p).astype("float32") if p.exists() else None

    def set_ndvi(self, cid, arr):
        np.save(self.dir / f"n_{_slug(cid)}.npy", arr.astype("float16"))

    def get_rgb(self, cid):
        p = self.dir / f"rgb_{_slug(cid)}.npy"
        return np.load(p) if p.exists() else None

    def set_rgb(self, cid, arr):
        np.save(self.dir / f"rgb_{_slug(cid)}.npy", arr)


# =============================================================================
# Persistência (desvio relativo à mediana do talhão, por janela)
# =============================================================================
def suavizar(a, mask, sigma=1.0):
    m = (mask & np.isfinite(a)).astype("float32")
    num = ndi.gaussian_filter(np.where(m > 0, a, 0).astype("float32"), sigma, mode="nearest")
    den = ndi.gaussian_filter(m, sigma, mode="nearest")
    out = np.where(den > 0.05, num / np.maximum(den, 1e-6), np.nan).astype("float32")
    out[~mask] = np.nan
    return out


def persistencia(P, mask, p: Params):
    total = max(int(mask.sum()), 1)
    D, usadas = [], []
    for k in range(P.shape[0]):
        v = P[k]
        ok = mask & np.isfinite(v)
        if ok.sum() / total < p.cobertura_min:
            continue
        b = float(np.median(v[ok]) if p.base == "mediana" else np.mean(v[ok]))
        if b < p.ndvi_min_janela:
            continue
        d = np.full(v.shape, np.nan, dtype="float32")
        d[ok] = (v[ok] - b) / b
        D.append(suavizar(d, mask, 1.0))          # suaviza o desvio ANTES de aplicar limiares
        usadas.append(k)
    if len(usadas) < 4:
        return None
    D = np.stack(D)
    n = np.isfinite(D).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mu, sd = np.nanmean(D, axis=0), np.nanstd(D, axis=0)
        nr = max(1, min(p.n_recente, len(usadas) // 2))
        rec, ant = np.nanmean(D[-nr:], axis=0), np.nanmean(D[:-nr], axis=0)
    lim = p.limiar_anomalia
    ult = D[-nr:]
    recneg = np.sum(ult < -lim, axis=0) / np.maximum(np.isfinite(ult).sum(axis=0), 1)   # fração das últimas janelas em queda
    den = np.maximum(n, 1)
    fneg = (np.sum(D < -lim, axis=0) / den).astype("float32")
    fpos = (np.sum(D > lim, axis=0) / den).astype("float32")
    inval = n < max(3, math.ceil(0.5 * len(usadas)))
    feats = dict(mu=mu, sigma=sd, fneg=fneg, fpos=fpos, rec=rec, ant=ant, recneg=recneg)
    for k, a in feats.items():
        a = a.astype("float32")
        a[inval] = np.nan
        feats[k] = a
    return dict(usadas=usadas, D=D, n=n, **feats)


def limpar_mascara(m, min_px):
    if not m.any():
        return m
    est = np.ones((3, 3), bool)
    m = ndi.binary_closing(ndi.binary_opening(m, est), est)
    m = ndi.gaussian_filter(m.astype("float32"), 1.0) > 0.5          # arredonda cantos de pixel
    lab, nc = ndi.label(m, structure=est)
    if nc:
        tam = ndi.sum(m, lab, index=np.arange(1, nc + 1))
        m = np.isin(lab, np.where(tam >= min_px)[0] + 1)
    bur, nb = ndi.label(~m)
    if nb:
        tam = ndi.sum(~m, bur, index=np.arange(1, nb + 1))
        m = m | np.isin(bur, np.where(tam < min_px)[0] + 1)
    return m


def camadas_estabilidade_interesse(f, mask, p: Params, px_ha):
    lim = p.limiar_anomalia
    min_px = max(1, int(round(p.area_min_camada_ha / px_ha)))
    ok = np.isfinite(f["mu"]) & mask
    persist = ok & (f["fneg"] >= p.f_interesse)
    # queda recente: abaixo do limiar em >= 2/3 das últimas janelas, e ainda normal nas anteriores
    recente = ok & np.isfinite(f["recneg"]) & np.isfinite(f["ant"]) & (f["recneg"] >= 0.66) & (f["ant"] > -lim / 2) & ~persist
    estavel = ok & (np.abs(f["mu"]) <= lim / 2) & ((f["fneg"] + f["fpos"]) <= p.f_estavel) & ~persist & ~recente
    persist, recente, estavel = (limpar_mascara(m, min_px) & mask for m in (persist, recente, estavel))
    recente &= ~persist
    estavel &= ~persist & ~recente
    return dict(estavel=estavel, persistente=persist, recente=recente)


# =============================================================================
# Zoneamento com restrições flexíveis de área
# =============================================================================
def _adjacencias(seg2d):
    pares = set()
    for a, b in ((seg2d[:, :-1], seg2d[:, 1:]), (seg2d[:-1, :], seg2d[1:, :])):
        m = (a > 0) & (b > 0) & (a != b)
        pares |= set(zip(np.minimum(a[m], b[m]).tolist(), np.maximum(a[m], b[m]).tolist()))
    return pares


def _arestas(seg2d, K):
    """Arestas compartilhadas entre segmentos e perímetro de cada segmento (em arestas de pixel)."""
    sh = {}
    for a, b in ((seg2d[:, :-1], seg2d[:, 1:]), (seg2d[:-1, :], seg2d[1:, :])):
        m = (a > 0) & (b > 0) & (a != b)
        for i, j in zip(np.minimum(a[m], b[m]).tolist(), np.maximum(a[m], b[m]).tolist()):
            sh[(i, j)] = sh.get((i, j), 0) + 1
    per = np.zeros(K + 1)
    ny, nx = seg2d.shape
    pad = np.pad(seg2d, 1)
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        nb = pad[1 + dy:1 + dy + ny, 1 + dx:1 + dx + nx]
        m = (seg2d > 0) & (nb != seg2d)
        np.add.at(per, seg2d[m], 1)
    return sh, per


def fundir_rag(seg2d, seg, Xc, N, px_ha, amin, amax, lam, ss_tot, piso_ha, lam_comp=0.0):
    """Fusão gananciosa de segmentos contíguos: Ward + penalidade flexível de área
    + penalidade de compacidade (índice de Polsby-Popper em pixels)."""
    K = int(seg.max())
    d = Xc.shape[1]
    n = np.bincount(seg, minlength=K + 1).astype(float)
    S = np.zeros((K + 1, d))
    np.add.at(S, seg, Xc)
    area = n * px_ha
    sh, per = _arestas(seg2d, K)
    adj = [set() for _ in range(K + 1)]
    for a, b in sh:
        adj[a].add(b)
        adj[b].add(a)
    pmax = lambda a: (max(0.0, a - amax) / amax) ** 2
    pmin = lambda a: (max(0.0, amin - a) / amin) ** 2
    pen = lambda a: pmax(a) + pmin(a)

    def penc(npx, p):                      # 1 = quadrado perfeito; penaliza abaixo de 0,6
        return max(0.0, 0.6 - (4 * math.pi * npx / (p * p)) / (math.pi / 4 * 4)) ** 2 if p > 0 else 0.0

    key = lambda a, b: (a, b) if a < b else (b, a)

    def custo(i, j, com_pen=True):
        dss = n[i] * n[j] / (n[i] + n[j]) * float(np.sum((S[i] / n[i] - S[j] / n[j]) ** 2))
        c = dss / ss_tot
        if com_pen:
            c += lam * (pen(area[i] + area[j]) - pen(area[i]) - pen(area[j]))
            if lam_comp:
                pn = per[i] + per[j] - 2 * sh.get(key(i, j), 0)
                c += lam_comp * (penc(n[i] + n[j], pn) - penc(n[i], per[i]) - penc(n[j], per[j]))
        return c

    vivos = set(range(1, K + 1))
    ver = [0] * (K + 1)
    pai = list(range(K + 1))
    heap = [(custo(i, j), 0, 0, i, j) for i in vivos for j in adj[i] if i < j]
    heapq.heapify(heap)

    def fundir(i, j):
        per[i] = per[i] + per[j] - 2 * sh.pop(key(i, j), 0)
        S[i] += S[j]; n[i] += n[j]; area[i] += area[j]
        vivos.discard(j); ver[i] += 1; pai[j] = i
        for k in adj[j]:
            if k != i:
                sh[key(i, k)] = sh.get(key(i, k), 0) + sh.pop(key(j, k), 0)
                adj[k].discard(j); adj[k].add(i); adj[i].add(k)
        adj[i].discard(j); adj[j] = set()
        for k in adj[i]:
            heapq.heappush(heap, (custo(i, k), ver[i], ver[k], i, k))

    while len(vivos) > N and heap:
        c, vi, vj, i, j = heapq.heappop(heap)
        if i in vivos and j in vivos and ver[i] == vi and ver[j] == vj:
            fundir(i, j)
    # piso duro: zonas inviáveis de amostrar são fundidas ao vizinho mais parecido
    while len(vivos) > 1:
        peq = min(vivos, key=lambda k: area[k])
        if area[peq] >= piso_ha or not adj[peq]:
            break
        viz = min(adj[peq], key=lambda k: custo(peq, k, com_pen=False))
        fundir(viz, peq) if viz != peq else None

    def raiz(i):
        while pai[i] != i:
            pai[i] = pai[pai[i]]
            i = pai[i]
        return i
    mapa = {r: idx + 1 for idx, r in enumerate(sorted(vivos))}
    return np.array([mapa[raiz(int(s))] for s in seg])


def _alocar(areas, N_total):
    areas = np.asarray(areas, float)
    N_total = max(int(N_total), len(areas))
    ideal = N_total * areas / areas.sum()
    aloc = np.maximum(1, np.floor(ideal)).astype(int)
    while aloc.sum() < N_total:
        aloc[np.argmax(ideal - aloc)] += 1
    while aloc.sum() > N_total and (aloc > 1).any():
        cand = np.where(aloc > 1)[0]
        aloc[cand[np.argmin((ideal - aloc)[cand])]] -= 1
    return aloc


def xy_padronizado(mask):
    Y, X = np.mgrid[0:mask.shape[0], 0:mask.shape[1]].astype("float32")
    out = []
    for a in (X, Y):
        v = a[mask]
        out.append(np.where(mask, (a - v.mean()) / (v.std() + 1e-9), 0))
    return np.stack(out, -1).astype("float32")


def zonear(X2d, XY2d, mask, px_ha, p: Params, ss_tot, n_total=None):
    """X2d: (ny,nx,d) covariáveis padronizadas e ponderadas; XY2d: (ny,nx,2) coordenadas
    padronizadas. Etapa 1: fragmentos compactos (Ward com peso espacial alto, contíguos).
    Etapa 2: fusão por semelhança com penalidades flexíveis de área e de compacidade.
    Devolve rótulos 1..N (0 fora da máscara)."""
    c = float(np.clip(p.compacidade, 0, 1))
    w_semente, w_fusao, lam_comp = 0.5 + 3.0 * c, 0.6 * c, 0.1 * c
    comp, nc = ndi.label(mask)           # 4-conectividade
    if nc == 0:
        raise RuntimeError("Máscara vazia.")
    areas = [float((comp == k).sum() * px_ha) for k in range(1, nc + 1)]
    N_total = int(n_total) if n_total else max(1, int(round(sum(areas) / p.densidade_ha)))
    aloc = _alocar(areas, N_total)
    out = np.zeros(mask.shape, int)
    base = 0
    for k in range(1, nc + 1):
        sub = comp == k
        Nc = int(aloc[k - 1]) if areas[k - 1] >= p.area_min_ha else 1
        ys, xs = np.where(sub)
        y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
        m = sub[y0:y1, x0:x1]
        F = X2d[y0:y1, x0:x1][m]
        XYc = XY2d[y0:y1, x0:x1][m]
        n = F.shape[0]
        if Nc <= 1 or n < 30:
            lab = np.ones(n, int)
        else:
            K0 = max(Nc, int(min(max(6 * Nc, Nc + 4), n // 6)))
            conn = grid_to_graph(m.shape[1], m.shape[0], mask=m)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                seg = AgglomerativeClustering(n_clusters=K0, connectivity=conn, linkage="ward").fit_predict(
                    np.hstack([F, XYc * w_semente])) + 1
            seg2d = np.zeros(m.shape, int)
            seg2d[m] = seg
            lab = fundir_rag(seg2d, seg, np.hstack([F, XYc * w_fusao]), Nc, px_ha, p.area_min_ha, p.area_max_ha,
                             p.rigidez_area, ss_tot, p.piso_ha, lam_comp)
        tmp = np.zeros(m.shape, int)
        tmp[m] = lab + base
        out[y0:y1, x0:x1] += tmp
        base += int(lab.max())
    return out


def preencher_vizinho(lab, alvo):
    if not (lab > 0).any():
        return lab
    idx = ndi.distance_transform_edt(lab == 0, return_distances=False, return_indices=True)
    cheio = lab[tuple(idx)]
    return np.where(alvo, cheio, 0)


def suavizar_rotulos(lab):
    ks = np.unique(lab[lab > 0])
    if len(ks) < 2:
        return lab
    memb = np.stack([ndi.gaussian_filter((lab == k).astype("float32"), 1.5) for k in ks])
    novo = ks[np.argmax(memb, axis=0)]
    novo[lab == 0] = 0
    for k in ks:
        cc, nc = ndi.label(novo == k)
        if nc > 1:
            tam = ndi.sum(novo == k, cc, index=np.arange(1, nc + 1))
            novo[(novo == k) & (cc != (np.argmax(tam) + 1))] = 0
    novo = preencher_vizinho(novo, lab > 0)
    return novo if len(np.unique(novo[novo > 0])) == len(ks) else lab


# =============================================================================
# Suavização dos contornos (sem frestas nem sobreposição)
# =============================================================================
def _sigma_suav(s):
    """Raio da suavização gaussiana, em pixels de trabalho."""
    return 0.8 + 1.2 * s


def _tol_suav(s, res):
    """Tolerância de simplificação (m) aplicada antes da suavização vetorial."""
    return 0.25 * res


def _sigma_vet(s, res):
    """Raio (m) da suavização vetorial das fronteiras, medido ao longo da curva."""
    return (0.5 + 1.1 * s) * res


def _fator_fino(grade, alvo_px=3_000_000):
    f = 4
    while f > 1 and (grade.ny * f) * (grade.nx * f) > alvo_px:
        f -= 1
    return f


def affine_fino(grade, f):
    return Affine(grade.res / f, 0, grade.x0, 0, -grade.res / f, grade.y0)


def particao_suave(lab, s, grade):
    """Suaviza a PARTIÇÃO inteira: a pertinência gaussiana de cada zona é interpolada numa grade f vezes
    mais fina e cada ponto fica com a zona de maior pertinência. Como é um argmax, as zonas continuam
    se encaixando (sem frestas nem sobreposição) e as fronteiras compartilhadas são suavizadas uma só vez."""
    f = _fator_fino(grade)
    ks = np.unique(lab[lab > 0])
    ny, nx = lab.shape
    perto = np.kron(ndi.binary_dilation(lab > 0, iterations=3), np.ones((f, f), bool))
    for tent in range(4):                                   # se alguma zona sumir, reduz a suavização
        sigma = _sigma_suav(s) * (0.7 ** tent)
        best = np.full((ny * f, nx * f), -1.0, dtype="float32")
        out = np.zeros((ny * f, nx * f), dtype="int32")
        for k in ks:
            m = ndi.gaussian_filter((lab == k).astype("float32"), sigma, mode="nearest")
            mf = ndi.zoom(m, f, order=3, mode="nearest", grid_mode=True)       # spline cúbica: contornos sem quinas
            up = mf > best
            best[up] = mf[up]
            out[up] = k
        out[~perto] = 0
        if len(np.unique(out[out > 0])) == len(ks):
            break
    for k in ks:                                            # um único bloco por zona
        cc, nc = ndi.label(out == k)
        if nc > 1:
            tam = ndi.sum(out == k, cc, index=np.arange(1, nc + 1))
            out[(out == k) & (cc != np.argmax(tam) + 1)] = 0
    return preencher_vizinho(out, perto & (out >= 0)), f


def mascara_suave(m, s, grade):
    f = _fator_fino(grade)
    mm = ndi.gaussian_filter(m.astype("float32"), _sigma_suav(s), mode="nearest")
    return ndi.zoom(mm, f, order=3, mode="nearest", grid_mode=True) > 0.5, f


# =============================================================================
# Polígonos
# =============================================================================
def poligonizar(rotulos, affine):
    out = {}
    for g, v in rasterio.features.shapes(rotulos.astype("int32"), mask=rotulos > 0, transform=affine, connectivity=4):
        out.setdefault(int(v), []).append(shape(g))
    return {k: unary_union(v) for k, v in out.items()}


def _poligonos_de(g):
    if g.is_empty:
        return []
    return [g] if g.geom_type == "Polygon" else [x for x in getattr(g, "geoms", []) if x.geom_type == "Polygon"]


def limpar_geom(g, min_m2=500.0):
    """Mantém só partes poligonais com área mínima; devolve Polygon/MultiPolygon ou None."""
    partes = [q for q in _poligonos_de(g) if q.area >= min_m2]
    if not partes:
        return None
    g = unary_union([q.buffer(0) for q in partes])           # garante geometria válida
    return g if g.geom_type in ("Polygon", "MultiPolygon") else None


TAPER = 0.5          # faixa (em múltiplos de sigma) onde o arco volta ao desenho original junto aos nós


def _reamostrar(xy, passo):
    d = np.hypot(*np.diff(xy, axis=0).T)
    s = np.r_[0.0, np.cumsum(d)]
    n = max(int(s[-1] / passo), 2) + 1
    t = np.linspace(0, s[-1], n)
    return np.c_[np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])], t, s[-1]


def suavizar_arco(xy, sigma_m, passo=2.0, aberto=False):
    """Suaviza uma curva (média gaussiana ao longo do comprimento). Em arcos abertos os extremos
    ficam fixos (peso decai até zero nas pontas); anéis isolados são suavizados em laço.
    `aberto=True` força o tratamento de extremos fixos (laço preso a um nó compartilhado)."""
    fechado = np.allclose(xy[0], xy[-1]) and not aberto
    P, t, L = _reamostrar(xy, passo)
    sg = max(sigma_m / passo, 0.5)
    if fechado:
        Q = ndi.gaussian_filter1d(P[:-1], sg, axis=0, mode="wrap")
        return np.vstack([Q, Q[:1]])
    Q = ndi.gaussian_filter1d(P, sg, axis=0, mode="nearest")
    w = np.clip(np.minimum(t, L - t) / (TAPER * sigma_m), 0, 1)
    w = (w * w * (3 - 2 * w))[:, None]                   # transição suave: extremos intactos
    return P + (Q - P) * w


def suavizar_rede(zonas, uni, sigma_m, passo=2.0):
    """Suaviza as FRONTEIRAS COMPARTILHADAS entre zonas como arcos (entre os pontos onde 3+ zonas se
    encontram, que ficam fixos) e reconstrói os polígonos. A borda do talhão não é alterada."""
    from shapely.geometry import LineString
    from shapely.ops import linemerge, polygonize
    if sigma_m <= 0 or not zonas:
        return zonas
    zs = list(zonas)
    linhas = []
    for g in zs:
        for q in _poligonos_de(g):
            linhas += [LineString(r.coords) for r in [q.exterior, *q.interiors]]
    rede = unary_union(linhas)
    arcos = linemerge(rede) if rede.geom_type == "MultiLineString" else rede
    arcos = list(getattr(arcos, "geoms", [arcos]))
    borda = uni.boundary
    chave = lambda p: (round(float(p[0]), 2), round(float(p[1]), 2))
    from collections import Counter
    toques = Counter()
    for a in arcos:
        toques[chave(a.coords[0])] += 1
        toques[chave(a.coords[-1])] += 1
    dados = []                                               # (xy original, curva suavizada ou None)
    for a in arcos:
        xy = np.asarray(a.coords)
        amostra = [a.interpolate(t, normalized=True) for t in (0.1, 0.5, 0.9)]
        if len(xy) < 4 or all(pt.distance(borda) < 0.5 for pt in amostra):
            dados.append((xy, None))                         # borda do talhão: intacta
        else:
            laco_preso = np.allclose(xy[0], xy[-1]) and toques[chave(xy[0])] > 2   # laço preso a um nó
            dados.append((xy, suavizar_arco(xy, sigma_m, passo, aberto=laco_preso)))
    # nós onde 3+ zonas se encontram: acomodam-se ao centro das curvas que chegam neles (sem "dentes")
    fixos, chegando = set(), {}
    for xy, Q in dados:
        for idx in (0, -1):
            n = chave(xy[idx])
            if Q is None:
                fixos.add(n)
                continue
            s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(Q, axis=0).T))]
            alvo = min(sigma_m, s[-1] / 2)                   # ponto da curva a ~sigma do nó
            j = int(np.searchsorted(s, alvo if idx == 0 else s[-1] - alvo))
            chegando.setdefault(n, []).append(Q[min(j, len(Q) - 1)])
    desloc = {}
    for n, pts in chegando.items():
        if n in fixos or len(pts) < 3:
            continue
        desloc[n] = np.mean(pts, axis=0) - np.array(n)

    def montar(fator):
        """Arcos com os nós deslocados em `fator` do valor ideal; polígonos reconstruídos."""
        novos = []
        for xy, Q in dados:
            if Q is None:
                novos.append(LineString(xy))
                continue
            s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(Q, axis=0).T))]
            L = s[-1]
            for idx, dist in ((0, s), (-1, L - s)):
                d = desloc.get(chave(xy[idx]))
                if d is not None and fator > 0:
                    w = np.clip(dist / min(1.5 * sigma_m, L), 0, 1)       # some exatamente na outra ponta
                    w = (1 - w * w * (3 - 2 * w))[:, None]
                    Q = Q + w * (fator * d)
            novos.append(LineString(Q).simplify(0.8))        # curva lisa, poucos vértices
        faces = list(polygonize(unary_union(novos)))
        saida = [[] for _ in zs]
        for f in faces:
            pt = f.representative_point()
            dono = next((i for i, g in enumerate(zs) if g.contains(pt)), None)
            if dono is None:
                dono = int(np.argmax([g.intersection(f).area for g in zs]))
            saida[dono].append(f)
        return saida

    saida = None
    for fator in (1.0, 0.5, 0.25, 0.0):                      # se arcos se cruzarem, desloca menos
        tent = montar(fator)
        if not all(tent):
            continue
        areas = [sum(f.area for f in s) for s in tent]
        ok_tot = abs(sum(areas) - sum(g.area for g in zs)) < 0.005 * sum(g.area for g in zs)
        ok_zona = all(0.6 < a / g.area < 1.4 for a, g in zip(areas, zs))       # nenhuma zona distorcida
        if ok_tot and ok_zona:
            saida = tent
            break
    if saida is None:
        return zonas
    novas = [limpar_geom(unary_union(s)) for s in saida]
    if any(g is None for g in novas):
        return zonas
    return novas


def suavizar_poligono(g, sigma_m, passo=2.0):
    """Suaviza o(s) anel(éis) de um polígono isolado (camadas de estabilidade/interesse)."""
    from shapely.geometry import Polygon as _P
    partes = []
    for q in _poligonos_de(g):
        try:
            sh = _P(suavizar_arco(np.asarray(q.exterior.coords), sigma_m, passo),
                    [suavizar_arco(np.asarray(h.coords), sigma_m, passo) for h in q.interiors]).buffer(0)
        except Exception:
            sh = q
        sh = sh.simplify(0.8)
        partes.append(sh if (sh.is_valid and not sh.is_empty and sh.area > 0.5 * q.area) else q)
    return unary_union(partes)


# =============================================================================
# Resultado
# =============================================================================
@dataclass
class Resultado:
    nome: str
    epsg: int
    grade: Grade
    params: Params
    campo_ll: object
    zonas: list
    estavel: list
    interesse: list
    rgb: np.ndarray | None
    rgb_grade: Grade | None
    rgb_data: str | None
    rasters: dict
    trafego: Trafego
    janelas_usadas: list
    janelas_total: list
    avisos: list
    fonte: str

    def tabela(self):
        return pd.DataFrame([{k: v for k, v in z.items() if k not in ("geom_utm", "geom_ll")} for z in self.zonas])

    def resumo(self):
        a = [z["area_ha"] for z in self.zonas]
        tot = sum(a)
        return dict(n_zonas=len(a), area_total_ha=tot, ha_por_amostra=tot / max(len(a), 1),
                    area_min=min(a), area_max=max(a),
                    fora_dos_limites=sum(1 for x in a if x < self.params.area_min_ha or x > self.params.area_max_ha),
                    pct_estavel=100 * sum(g["area_ha"] for g in self.estavel) / tot,
                    pct_interesse=100 * sum(g["area_ha"] for g in self.interesse) / tot)


class _Ctx:
    pass


# =============================================================================
# Pipeline por talhão
# =============================================================================
def _adquirir(ctx, p, prog):
    """Séries de NDVI p90 por janela, lendo o mínimo necessário (triagem SCL -> NDVI)."""
    ini, fim = date.fromisoformat(p.inicio), date.fromisoformat(p.fim)
    janelas = zp.agrupar_por_janela(ctx.cenas, ini, fim)
    P, rotulos = [], []
    pool_n = p.workers

    def triar(c):
        v = ctx.cache.tri.get(c.id)
        if v is not None:
            with ctx.lock:
                ctx.traf.triagens_cache += 1
            return v
        v = ctx.fonte.claro(c, ctx.grade20, ctx.mask20)
        ctx.cache.set_tri(c.id, v)
        with ctx.lock:
            ctx.traf.triagens_rede += 1
        return v

    def obter_ndvi(c):
        a = ctx.cache.get_ndvi(c.id)
        if a is not None:
            with ctx.lock:
                ctx.traf.ndvi_cache += 1
            return a
        a = ctx.fonte.ndvi(c, ctx.grade)
        ctx.cache.set_ndvi(c.id, a)
        with ctx.lock:
            ctx.traf.ndvi_rede += 1
        return a.astype("float16").astype("float32")

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(pool_n) as ex:
        for i, (rot, cs) in enumerate(janelas):
            prog(0.15 + 0.5 * i / max(len(janelas), 1), f"Janela {rot}: triagem de nuvem")
            cs = sorted(cs, key=lambda c: c.nuvem)
            boas = []
            for j in range(0, len(cs), pool_n):
                lote = cs[j:j + pool_n]
                for c, v in zip(lote, ex.map(triar, lote)):
                    if v >= p.claro_min:
                        boas.append((v, c))
                    else:
                        ctx.traf.cenas_descartadas_nuvem += 1
                if len(boas) >= p.obs_por_janela:
                    break
            boas = [c for _, c in sorted(boas, key=lambda t: -t[0])[:p.obs_por_janela]]
            if len(boas) < p.obs_min_janela:
                ctx.avisos.append(f"Janela {rot} descartada: só {len(boas)} cena(s) sem nuvem.")
                continue
            arrs = list(ex.map(obter_ndvi, boas))
            pj, _ = zp.p90_janela(np.stack(arrs), p.obs_min_pixel)
            P.append(pj)
            rotulos.append(rot)
    ctx.cache.salvar_tri()
    if not P:
        raise RuntimeError("Nenhuma janela com cenas suficientes.")
    return np.stack(P), rotulos


def _relevo_na_grade(ctx, p, uni_utm, mask_full):
    gdem = grade_alinhada(uni_utm.bounds, ctx.epsg, RES_DEM, 1500, alinhar=RES_DEM)
    z = ctx.fonte.dem(ctx.bbox_dem, gdem)
    ctx.traf.dem_rede += 1
    rel = zp.derivar_relevo(z, RES_DEM)
    rm = rasterio.features.rasterize([(uni_utm, 1)], out_shape=gdem.shape, transform=gdem.affine,
                                     fill=0, dtype="uint8", all_touched=True).astype(bool)
    rel["elev_rel"] = rel["z"] - float(np.mean(rel["z"][rm]))
    saida = {}
    for k in ("elev_rel", "decliv_pct", "twi", "tpi200"):
        dst = np.full(ctx.grade.shape, np.nan, dtype="float32")
        rasterio.warp.reproject(rel[k].astype("float32"), dst, src_transform=gdem.affine,
                                src_crs=f"EPSG:{ctx.epsg}", dst_transform=ctx.grade.affine,
                                dst_crs=f"EPSG:{ctx.epsg}", resampling=Resampling.bilinear,
                                src_nodata=np.nan, dst_nodata=np.nan)
        saida[k] = dst
    return saida


def _preencher_nan(a, mask):
    ok = mask & np.isfinite(a)
    if ok.all() or not ok.any():
        return a
    idx = ndi.distance_transform_edt(~ok, return_distances=False, return_indices=True)
    return np.where(mask, a[tuple(idx)], np.nan)


def _escolher_rgb(ctx, p):
    fim = date.fromisoformat(p.fim)
    recentes = [c for c in ctx.cenas if (fim - c.datetime.date()).days <= p.dias_rgb]
    if not recentes:
        recentes = sorted(ctx.cenas, key=lambda c: c.datetime, reverse=True)[:12]
    recentes = sorted(recentes, key=lambda c: c.nuvem)[:8]
    melhor, mv = None, -1.0
    for c in recentes:
        v = ctx.cache.tri.get(c.id)
        if v is None:
            v = ctx.fonte.claro(c, ctx.grade20, ctx.mask20)
            ctx.cache.set_tri(c.id, v)
            ctx.traf.triagens_rede += 1
        else:
            ctx.traf.triagens_cache += 1
        if v > mv:
            melhor, mv = c, v
        if v >= 0.95:
            break
    ctx.cache.salvar_tri()
    return melhor


def analisar_talhao(nome, geom_ll, fonte, p: Params, cenas=None, prog=lambda f, m: None):
    t0 = time.time()
    ctx = _Ctx()
    ctx.fonte, ctx.traf, ctx.avisos = fonte, Trafego(), []
    ctx.lock = threading.Lock()
    ini, fim = p.inicio, p.fim
    bbox = tuple(geom_ll.buffer(0.001).bounds)
    ctx.bbox_dem = tuple(geom_ll.buffer(0.02).bounds)
    prog(0.02, "Buscando cenas (só metadados)")
    if cenas is None:
        cenas = fonte.buscar(bbox, ini, fim, p.nuvem_max_cena)
    if not cenas:
        raise RuntimeError("Nenhuma cena Sentinel-2 encontrada para o talhão.")
    epsg = max(set(c.epsg for c in cenas), key=[c.epsg for c in cenas].count)
    ctx.cenas = [c for c in cenas if c.epsg == epsg]
    ctx.epsg = epsg
    ctx.traf.cenas_encontradas = len(ctx.cenas)

    f = zp.func_transf("EPSG:4326", f"EPSG:{epsg}")
    uni = shp_transform(f, geom_ll)
    interior = uni.buffer(-p.buffer_borda_m)
    if interior.is_empty:
        ctx.avisos.append("Talhão pequeno: recuo de borda desativado.")
        interior = uni
    ctx.grade = grade_alinhada(uni.bounds, epsg, p.res_m, 40, alinhar=20)
    ctx.grade20 = grade_alinhada(uni.bounds, epsg, 20, 40, alinhar=20)
    rast = lambda g, gr, at=False: rasterio.features.rasterize([(g, 1)], out_shape=gr.shape, transform=gr.affine,
                                                              fill=0, dtype="uint8", all_touched=at).astype(bool)
    mask_full = rast(uni, ctx.grade, True)
    mask_int = rast(interior, ctx.grade)
    ctx.mask20 = rast(interior, ctx.grade20)
    chave = hashlib.md5(f"{geom_ll.wkb_hex}|{p.res_m}|{p.buffer_borda_m}|{epsg}|{getattr(fonte, 'colecao', 'syn')}".encode()).hexdigest()[:12]
    ctx.cache = Cache(p.cache_dir, f"{_slug(nome)}_{chave}")

    with MedidorRede() as med:
        P, rot = _adquirir(ctx, p, prog)
        prog(0.7, "Relevo (Copernicus DEM)")
        relevo = _relevo_na_grade(ctx, p, uni, mask_full)
        prog(0.78, "Imagem RGB (uma cena)")
        rgb, rgb_data, rgb_grade = None, None, None
        try:
            c_rgb = _escolher_rgb(ctx, p)
            rgb_grade = grade_alinhada(uni.bounds, epsg, 10, 60, alinhar=20)
            rgb = ctx.cache.get_rgb(c_rgb.id)
            if rgb is not None:
                ctx.traf.rgb_cache += 1
            else:
                rgb = fonte.rgb(c_rgb, rgb_grade)
                ctx.cache.set_rgb(c_rgb.id, rgb)
                ctx.traf.rgb_rede += 1
            rgb_data = c_rgb.datetime.date().isoformat()
        except Exception as e:  # RGB é opcional: não derruba a análise
            ctx.avisos.append(f"RGB indisponível: {e!r}")
    ctx.traf.bytes_medidos = med.bytes

    prog(0.85, "Zoneamento")
    met = persistencia(P, mask_int, p)
    if met is None:
        raise RuntimeError("Menos de 4 janelas com lavoura e cobertura suficientes; "
                           "reduza 'NDVI mínimo da janela' ou amplie o período.")
    px_ha = (p.res_m ** 2) / 1e4
    camp = {k: met[k] for k in ("mu", "sigma", "fneg", "fpos")}
    rel = {k: suavizar(relevo[k], mask_int) for k in relevo}
    cols = [("mu", camp["mu"], 1.0), ("sigma", camp["sigma"], 0.5), ("fneg", camp["fneg"], 0.8),
            ("fpos", camp["fpos"], 0.5), ("elev_rel", rel["elev_rel"], p.peso_relevo),
            ("decliv_pct", rel["decliv_pct"], 0.6 * p.peso_relevo), ("twi", rel["twi"], 0.8 * p.peso_relevo),
            ("tpi200", rel["tpi200"], 0.6 * p.peso_relevo)]
    camadas = []
    for nm, a, w in cols:
        a = _preencher_nan(a, mask_int)
        if not np.isfinite(a[mask_int]).all() or w <= 0:
            continue
        v = a[mask_int]
        sd = v.std()
        if sd < 1e-9:
            continue
        camadas.append(np.where(mask_int, (a - v.mean()) / sd * w, 0).astype("float32"))
    X2d = np.stack(camadas, axis=-1)
    Xv = X2d[mask_int]
    ss_tot = float(np.sum((Xv - Xv.mean(0)) ** 2)) + 1e-9
    n_alvo = max(1, int(round(zp.area_geodesica_ha(geom_ll) / p.densidade_ha)))   # sobre a área TOTAL do talhão
    lab_int = zonear(X2d, xy_padronizado(mask_int), mask_int, px_ha, p, ss_tot, n_alvo)
    prog(0.93, "Polígonos")
    if p.suavizacao > 0:
        lab_full = preencher_vizinho(lab_int, mask_full)
        lab_fino, fator = particao_suave(lab_full, p.suavizacao, ctx.grade)
        polis = poligonizar(lab_fino, affine_fino(ctx.grade, fator))
        tol_zona = _tol_suav(p.suavizacao, ctx.grade.res)
    else:
        lab_full = preencher_vizinho(suavizar_rotulos(lab_int), mask_full)
        polis = poligonizar(lab_full, ctx.grade.affine)
        tol_zona = 1.0 * ctx.grade.res
    ids = sorted(polis)
    if hasattr(shapely, "coverage_simplify"):          # Shapely >= 2.1: preserva as fronteiras entre zonas
        simp = shapely.coverage_simplify([polis[i] for i in ids], tol_zona)
    else:                                               # fallback: pode deixar frestas mínimas entre zonas
        ctx.avisos.append("Shapely < 2.1: zonas podem ter frestas/geometria inválida. Atualize: pip install -U 'shapely>=2.1'")
        simp, ocupado = [], None
        for i in ids:
            g = polis[i].simplify(tol_zona, preserve_topology=True)
            g = g if ocupado is None else g.difference(ocupado)       # elimina sobreposições
            simp.append(g)
            ocupado = g if ocupado is None else ocupado.union(g)
    zcov = list(simp)                         # fronteiras compartilhadas idênticas (cobertura)
    if p.suavizacao > 0:
        if hasattr(shapely, "coverage_simplify"):
            bruto, externo = zcov, unary_union(zcov)
            sig = _sigma_vet(p.suavizacao, ctx.grade.res)
            for k, fr in enumerate((1.0, 0.6, 0.35)):               # zonas estreitas: tenta com raio menor
                zcov = suavizar_rede(bruto, externo, fr * sig)
                if zcov is not bruto:
                    if k:
                        ctx.avisos.append(f"Suavização das fronteiras reduzida a {fr:.0%} (há zonas estreitas).")
                    break
            else:
                ctx.avisos.append("Suavização das fronteiras entre zonas não aplicada (zonas muito estreitas); contornos mantidos como gerados.")
        else:
            ctx.avisos.append("Suavização das fronteiras exige Shapely >= 2.1 (pip install -U 'shapely>=2.1'); contornos mantidos como gerados.")
    zonas_utm = []
    for g in zcov:
        g = limpar_geom(g.intersection(uni))                          # só agora recorta pelo limite real do talhão
        if g is not None:
            zonas_utm.append(g)
    # numeração estável: norte->sul em faixas, depois oeste->leste
    faixa = 3 * ctx.grade.res * 8
    if p.suavizacao > 0:
        bruto = zonas_utm
        sig = _sigma_vet(p.suavizacao, ctx.grade.res)
        for k, fr in enumerate((1.0, 0.6, 0.35)):               # zonas estreitas: tenta com raio menor
            zonas_utm = suavizar_rede(bruto, uni, fr * sig)
            if zonas_utm is not bruto:
                if k:
                    ctx.avisos.append(f"Suavização das fronteiras reduzida a {fr:.0%} (há zonas estreitas).")
                break
        else:
            ctx.avisos.append("Suavização das fronteiras entre zonas não aplicada (zonas muito estreitas); contornos mantidos como gerados.")
    zonas_utm.sort(key=lambda g: (-round(g.representative_point().y / faixa), g.representative_point().x))
    f_ll = zp.func_transf(f"EPSG:{epsg}", "EPSG:4326")
    cam = camadas_estabilidade_interesse(
        dict(mu=camp["mu"], fneg=camp["fneg"], fpos=camp["fpos"], recneg=met["recneg"], ant=met["ant"]), mask_int, p, px_ha)

    def vetor(mask):
        if not mask.any():
            return []
        if p.suavizacao > 0:
            mf, fator = mascara_suave(mask, p.suavizacao, ctx.grade)
            g = unary_union(list(poligonizar(mf.astype(int), affine_fino(ctx.grade, fator)).values()))
            g = shapely.simplify(g, tol_zona)
        else:
            g = unary_union(list(poligonizar(mask.astype(int), ctx.grade.affine).values()))
            g = shapely.simplify(g, 0.6 * ctx.grade.res)
        if p.suavizacao > 0:
            g = suavizar_poligono(g, _sigma_vet(p.suavizacao, ctx.grade.res))
        g = limpar_geom(g.intersection(uni), 0.5 * p.area_min_camada_ha * 1e4)
        return _poligonos_de(g) if g is not None else []

    def com_attr(geoms_utm, **extra):
        out = []
        for g in geoms_utm:
            gl = shp_transform(f_ll, g)
            out.append(dict(geom_utm=g, geom_ll=gl, area_ha=zp.area_geodesica_ha(gl), **extra))
        return out

    interesse = com_attr(vetor(cam["persistente"]), tipo="queda persistente") + \
        com_attr(vetor(cam["recente"]), tipo="queda recente")
    est_utm = vetor(cam["estavel"])
    if interesse:                       # as camadas nunca se sobrepõem
        g_int = unary_union([g["geom_utm"] for g in interesse])
        est_utm = [q for g in est_utm for q in _poligonos_de(limpar_geom(g.difference(g_int), 0.5 * p.area_min_camada_ha * 1e4) or Polygon())]
    estavel = com_attr(est_utm, tipo="estavel")

    zonas = []
    lab_pix = {}
    for n_, g in enumerate(zonas_utm, 1):
        gm = rast(g, ctx.grade) & mask_int
        gl = shp_transform(f_ll, g)
        st = {}
        for k, a in dict(mu=camp["mu"], fneg=camp["fneg"], fpos=camp["fpos"], elev_rel=rel["elev_rel"],
                         decliv_pct=rel["decliv_pct"], twi=rel["twi"]).items():
            v = a[gm]
            v = v[np.isfinite(v)]
            st[k] = float(v.mean()) if v.size else np.nan
        tot = max(int(gm.sum()), 1)
        zonas.append(dict(zona=n_, nome=f"{nome}_Z{n_:02d}", area_ha=zp.area_geodesica_ha(gl),
                          lon=gl.representative_point().x, lat=gl.representative_point().y,
                          pct_estavel=100 * float((gm & cam["estavel"]).sum()) / tot,
                          pct_interesse=100 * float((gm & (cam["persistente"] | cam["recente"])).sum()) / tot,
                          **st, geom_utm=g, geom_ll=gl))
        lab_pix[n_] = gm
    rot_final = np.zeros(ctx.grade.shape, int)
    for n_, gm in lab_pix.items():
        rot_final[gm] = n_

    usadas = [rot[k] for k in met["usadas"]]
    prog(1.0, "Pronto")
    return Resultado(nome=nome, epsg=epsg, grade=ctx.grade, params=p, campo_ll=geom_ll, zonas=zonas,
                     estavel=estavel, interesse=interesse, rgb=rgb, rgb_grade=rgb_grade, rgb_data=rgb_data,
                     rasters=dict(mu=camp["mu"], fneg=camp["fneg"], fpos=camp["fpos"], rec=met["rec"], ant=met["ant"], recneg=met["recneg"], rotulos=rot_final,
                                  mask_int=mask_int, campo_utm=uni, **{f"c_{k}": v for k, v in cam.items()}),
                     trafego=ctx.traf, janelas_usadas=usadas, janelas_total=rot, avisos=ctx.avisos,
                     fonte=getattr(fonte, "nome", "?"))


def analisar_fazenda(talhoes, fonte, p: Params, prog=lambda f, m: None, erros=None):
    """talhoes: lista de (nome, geom_ll). Uma única busca STAC para a fazenda inteira.
    Se `erros` (lista) for passada, falhas por talhão são registradas e a análise continua."""
    allb = unary_union([g for _, g in talhoes]).buffer(0.001).bounds
    cenas_f = None
    if (allb[2] - allb[0]) < 0.4 and (allb[3] - allb[1]) < 0.4:
        prog(0.0, "Busca STAC única para a fazenda")
        cenas_f = fonte.buscar(allb, p.inicio, p.fim, p.nuvem_max_cena)
    out = []
    for i, (nome, g) in enumerate(talhoes):
        cenas = None
        if cenas_f is not None:
            minx, miny, maxx, maxy = g.buffer(0.001).bounds
            cenas = [c for c in cenas_f if c.item is None or _cena_intersecta(c, (minx, miny, maxx, maxy))]
        try:
            out.append(analisar_talhao(nome, g, fonte, p, cenas, lambda f, m, i=i, nome=nome: prog((i + f) / len(talhoes), f"{nome}: {m}")))
        except Exception as e:
            if erros is None:
                raise
            erros.append((nome, str(e)))
    return out


def talhao_exemplo():
    """Pivô de ~120 ha (para demonstração)."""
    c = zp.crs_local(shapely.Point(-47.9, -14.7))
    g = shp_transform(zp.func_transf(c, "EPSG:4326"), shapely.Point(0, 0).buffer(620, 64))
    return [("Pivô exemplo", g)]


def _cena_intersecta(c, bbox):
    b = c.item.bbox
    return not (b[2] < bbox[0] or b[0] > bbox[2] or b[3] < bbox[1] or b[1] > bbox[3])


# =============================================================================
# Entrada de dados
# =============================================================================
_NOMES = ("name", "nome", "talhao", "talhão", "field_name", "field", "id")


def ler_talhoes(caminhos, unir_por_arquivo=False):
    """Lê KML/KMZ/shapefile/GeoJSON. Devolve [(nome, geometria WGS84)]."""
    talhoes = []
    for c in caminhos:
        c = Path(c)
        if c.suffix.lower() == ".kmz":
            d = tempfile.mkdtemp(prefix="kmz_")
            extrair_zip_seguro(c, d)
            kml = next(Path(d).rglob("*.kml"))
            c = kml
        try:
            g = gpd.read_file(c)
        except Exception:
            g = gpd.read_file(c, driver="KML")
        if g.crs is None:
            g = g.set_crs(4326)
        g = g.to_crs(4326)
        g = g[g.geometry.notna() & g.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
        g["geometry"] = shapely.force_2d(g.geometry.values)
        g["geometry"] = g.geometry.make_valid()
        if g.empty:
            continue
        if unir_por_arquivo:
            talhoes.append((c.stem, unary_union(list(g.geometry))))
            continue
        col = next((k for k in g.columns if k.lower() in _NOMES), None)
        for i, (_, r) in enumerate(g.iterrows(), 1):
            nome = str(r[col]).strip() if col and str(r[col]).strip() not in ("", "None", "nan") else f"{c.stem}_{i:02d}"
            talhoes.append((nome, r.geometry))
    vistos, final = {}, []
    for n, g in talhoes:
        vistos[n] = vistos.get(n, 0) + 1
        final.append((n if vistos[n] == 1 else f"{n}_{vistos[n]}", g))
    return final


# =============================================================================
# KML
# =============================================================================
_PALETA = ["e6194b", "3cb44b", "ffe119", "4363d8", "f58231", "911eb4", "46f0f0", "f032e6", "bcf60c", "fabebe",
           "008080", "e6beff", "9a6324", "fffac8", "800000", "aaffc3", "808000", "ffd8b1", "000075", "808080"]


def _kml_cor(hex_rgb, alpha):
    r, g, b = hex_rgb[0:2], hex_rgb[2:4], hex_rgb[4:6]
    return f"{alpha}{b}{g}{r}"                       # KML usa aabbggrr


def _coords(ring):
    return " ".join(f"{x:.7f},{y:.7f},0" for x, y in ring.coords)


def _kml_poligono(g):
    out = []
    for q in _poligonos_de(g):
        q = orient(q, 1.0)
        s = f"<Polygon><outerBoundaryIs><LinearRing><coordinates>{_coords(q.exterior)}</coordinates></LinearRing></outerBoundaryIs>"
        for h in q.interiors:
            s += f"<innerBoundaryIs><LinearRing><coordinates>{_coords(h)}</coordinates></LinearRing></innerBoundaryIs>"
        out.append(s + "</Polygon>")
    return out[0] if len(out) == 1 else f"<MultiGeometry>{''.join(out)}</MultiGeometry>"


def _kml_doc(nome, corpo, estilos=""):
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
            f"<name>{escape(nome)}</name>{estilos}{corpo}</Document></kml>")


def kml_zonas(resultados, nome_doc="Zonas de amostragem direcionada", limpo=True, padrao_nome="{talhao}_Z{n:02d}"):
    """KML para upload: um Placemark (Polygon) por zona. `limpo` = só nome e geometria."""
    corpo, estilos = "", ""
    if not limpo:
        for i, c in enumerate(_PALETA):
            estilos += (f'<Style id="z{i}"><LineStyle><color>{_kml_cor("ffffff", "ff")}</color><width>2</width></LineStyle>'
                        f'<PolyStyle><color>{_kml_cor(c, "88")}</color></PolyStyle></Style>')
    for r in resultados:
        corpo += f"<Folder><name>{escape(r.nome)}</name>"
        for z in r.zonas:
            corpo += f"<Placemark><name>{escape(padrao_nome.format(talhao=r.nome, n=z['zona']))}</name>"
            if not limpo:
                desc = "Área %.1f ha | estável %.0f%% | interesse %.0f%%" % (z["area_ha"], z["pct_estavel"], z["pct_interesse"])
                corpo += ("<styleUrl>#z%d</styleUrl><description>%s</description>"
                          % ((z["zona"] - 1) % len(_PALETA), escape(desc)))
            corpo += _kml_poligono(z["geom_ll"]) + "</Placemark>"
        corpo += "</Folder>"
    return _kml_doc(nome_doc, corpo, estilos)


def kml_camadas(resultados, nome_doc="Camadas de análise (estabilidade e interesse)"):
    estilos = "".join(
        f'<Style id="{i}"><LineStyle><color>{_kml_cor(c, "ff")}</color><width>2</width></LineStyle>'
        f'<PolyStyle><color>{_kml_cor(c, "66")}</color></PolyStyle></Style>'
        for i, c in (("estavel", "2ca02c"), ("persistente", "d62728"), ("recente", "ff7f0e"), ("campo", "ffffff")))
    corpo = ""
    for r in resultados:
        corpo += f"<Folder><name>{escape(r.nome)}</name>"
        corpo += f"<Placemark><name>Limite</name><styleUrl>#campo</styleUrl>{_kml_poligono(r.campo_ll)}</Placemark>"
        for rot, itens, sid in (("Zonas de estabilidade", r.estavel, "estavel"),
                                ("Zonas de interesse — queda persistente", [g for g in r.interesse if g["tipo"] == "queda persistente"], "persistente"),
                                ("Zonas de interesse — queda recente", [g for g in r.interesse if g["tipo"] == "queda recente"], "recente")):
            corpo += f"<Folder><name>{rot}</name>"
            for k, g in enumerate(itens, 1):
                corpo += (f"<Placemark><name>{escape(rot.split('—')[-1].strip())} {k} ({g['area_ha']:.1f} ha)</name>"
                          f"<styleUrl>#{sid}</styleUrl>{_kml_poligono(g['geom_ll'])}</Placemark>")
            corpo += "</Folder>"
        corpo += "</Folder>"
    return _kml_doc(nome_doc, corpo, estilos)


def tabela_zonas(resultados):
    dfs = []
    for r in resultados:
        t = r.tabela()
        t.insert(0, "talhao", r.nome)
        dfs.append(t)
    return pd.concat(dfs, ignore_index=True)


# =============================================================================
# Visualização
# =============================================================================
def _desenhar(ax, g, **kw):
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath
    for q in _poligonos_de(g):
        verts, codes = [], []
        for ring in [q.exterior, *q.interiors]:
            xy = np.asarray(ring.coords)
            verts += list(xy)
            codes += [MPath.MOVETO] + [MPath.LINETO] * (len(xy) - 2) + [MPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MPath(verts, codes), **kw))


def painel_png(r: Resultado, dpi=110):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(12, 11))
    l, b, rr, t = r.rgb_grade.bounds() if r.rgb is not None else r.grade.bounds()
    campo = r.rasters["campo_utm"]
    for ax in axs.ravel():
        if r.rgb is not None:
            ax.imshow(np.moveaxis(r.rgb, 0, -1), extent=(l, rr, b, t), interpolation="bilinear")
        else:
            ax.set_facecolor("#dddddd")
        minx, miny, maxx, maxy = campo.bounds
        pad = 0.04 * max(maxx - minx, maxy - miny)
        ax.set_xlim(minx - pad, maxx + pad); ax.set_ylim(miny - pad, maxy + pad)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        _desenhar(ax, campo, fill=False, ec="white", lw=1.5)
    ax = axs[0, 0]
    ax.set_title(f"RGB Sentinel-2 ({r.rgb_data or 'sem imagem'})", fontsize=11)
    ax = axs[0, 1]
    for z in r.zonas:
        _desenhar(ax, z["geom_utm"], fc=_PALETA_MPL(z["zona"]), ec="white", lw=1.2, alpha=0.62)
        p_ = z["geom_utm"].representative_point()
        ax.text(p_.x, p_.y, f"Z{z['zona']:02d}\n{z['area_ha']:.1f} ha", ha="center", va="center", fontsize=7.5,
                color="black", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    ax.set_title(f"Zoneamento: {len(r.zonas)} zonas ({r.resumo()['ha_por_amostra']:.1f} ha/amostra)", fontsize=11)
    ax = axs[1, 0]
    for g in r.estavel:
        _desenhar(ax, g["geom_utm"], fc="#2ca02c", ec="#14501a", lw=1, alpha=0.6)
    ax.set_title(f"Estabilidade: {r.resumo()['pct_estavel']:.0f}% da área, sem desvios sobre a mediana", fontsize=11)
    ax = axs[1, 1]
    for g in r.interesse:
        cor = "#d62728" if g["tipo"] == "queda persistente" else "#ff7f0e"
        _desenhar(ax, g["geom_utm"], fc=cor, ec="black", lw=0.8, alpha=0.65)
    ax.set_title(f"Interesse (queda de vigor): {r.resumo()['pct_interesse']:.0f}% da área", fontsize=11)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="#d62728", label="queda persistente"), Patch(fc="#ff7f0e", label="queda recente")],
              loc="lower right", fontsize=8, framealpha=0.85)
    fig.suptitle(f"{r.nome} · {len(r.janelas_usadas)} janelas de safra · fonte: {r.fonte}", fontsize=12)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    plt.close(fig)
    return buf.getvalue()


def _PALETA_MPL(i):
    return "#" + _PALETA[(i - 1) % len(_PALETA)]


def rgb_wgs84(r: Resultado):
    g = r.rgb_grade
    l, b, rr, t = g.bounds()
    dst_t, w, h = rasterio.warp.calculate_default_transform(f"EPSG:{g.epsg}", "EPSG:4326", g.nx, g.ny, l, b, rr, t)
    out = np.zeros((4, h, w), "uint8")
    for i in range(3):
        rasterio.warp.reproject(r.rgb[i], out[i], src_transform=g.affine, src_crs=f"EPSG:{g.epsg}",
                                dst_transform=dst_t, dst_crs="EPSG:4326", resampling=Resampling.bilinear)
    rasterio.warp.reproject(np.full((g.ny, g.nx), 255, "uint8"), out[3], src_transform=g.affine,
                            src_crs=f"EPSG:{g.epsg}", dst_transform=dst_t, dst_crs="EPSG:4326",
                            resampling=Resampling.nearest)
    west, north = dst_t.c, dst_t.f
    return np.moveaxis(out, 0, -1), [[north + dst_t.e * h, west], [north, west + dst_t.a * w]]


def mapa_html(r: Resultado):
    """Mapa interativo (Leaflet/folium) sem mapa-base externo: o fundo é a própria RGB."""
    import folium
    c = r.campo_ll.centroid
    m = folium.Map(location=[c.y, c.x], zoom_start=15, tiles=None, control_scale=True, max_zoom=20)
    if r.rgb is not None:
        img, bounds = rgb_wgs84(r)
        folium.raster_layers.ImageOverlay(image=img, bounds=bounds, name=f"RGB ({r.rgb_data})", opacity=1.0, zindex=1).add_to(m)
    folium.GeoJson(mapping(r.campo_ll), name="Limite do talhão", style_function=lambda f: dict(color="white", weight=2, fillOpacity=0)).add_to(m)

    def camada(nome, feats, estilo, tip, show=True):
        fc = dict(type="FeatureCollection", features=[dict(type="Feature", geometry=mapping(g["geom_ll"]), properties=pr) for g, pr in feats])
        folium.GeoJson(fc, name=nome, style_function=estilo, show=show,
                       tooltip=folium.GeoJsonTooltip(fields=list(feats[0][1].keys()), aliases=tip) if feats else None).add_to(m)

    zf = [(z, dict(zona=f"Z{z['zona']:02d}", area=f"{z['area_ha']:.1f} ha", estavel=f"{z['pct_estavel']:.0f}%",
                   interesse=f"{z['pct_interesse']:.0f}%", cor=_PALETA_MPL(z["zona"]))) for z in r.zonas]
    camada("Clusterização (zonas propostas)", zf,
           lambda f: dict(color="white", weight=1.5, fillColor=f["properties"]["cor"], fillOpacity=0.5),
           ["Zona", "Área", "% estável", "% interesse", "cor"])
    ef = [(g, dict(tipo="estabilidade", area=f"{g['area_ha']:.1f} ha")) for g in r.estavel]
    if ef:
        camada("Zonas de estabilidade", ef, lambda f: dict(color="#14501a", weight=1, fillColor="#2ca02c", fillOpacity=0.55),
               ["Tipo", "Área"], show=False)
    inf = [(g, dict(tipo=g["tipo"], area=f"{g['area_ha']:.1f} ha")) for g in r.interesse]
    if inf:
        camada("Zonas de interesse (queda de vigor)", inf,
               lambda f: dict(color="black", weight=1, fillColor="#d62728" if f["properties"]["tipo"] == "queda persistente" else "#ff7f0e", fillOpacity=0.6),
               ["Tipo", "Área"], show=False)
    if getattr(r, "fonte", "").startswith("SINT"):
        m.get_root().html.add_child(folium.Element(
            '<div style="position:fixed;top:8px;left:55px;z-index:9999;background:#fff3cd;padding:4px 10px;'
            'border:1px solid #c9a227;font:600 13px sans-serif">DADOS SINTÉTICOS (demonstração)</div>'))
    folium.LayerControl(collapsed=False).add_to(m)
    m.fit_bounds([[r.campo_ll.bounds[1], r.campo_ll.bounds[0]], [r.campo_ll.bounds[3], r.campo_ll.bounds[2]]])
    return m.get_root().render()


# =============================================================================
# CLI
# =============================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description="Zoneamento direcionado TruBio (KML para upload).")
    ap.add_argument("entradas", nargs="*")
    ap.add_argument("--kml", default="zonas_propostas.kml")
    ap.add_argument("--camadas", default=None, help="KML extra com estabilidade/interesse")
    ap.add_argument("--painel", default=None, help="PNG com as 4 vistas (do primeiro talhão)")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--demo", action="store_true", help="dados sintéticos, sem rede")
    ap.add_argument("--unir-por-arquivo", action="store_true", help="um talhão por arquivo (polígonos unidos)")
    ap.add_argument("--selftest", action="store_true")
    for nome, tipo in (("densidade_ha", float), ("area_min_ha", float), ("area_max_ha", float), ("rigidez_area", float),
                       ("res_m", int), ("obs_por_janela", int), ("inicio", str), ("fim", str), ("cache_dir", str)):
        ap.add_argument(f"--{nome.replace('_', '-')}", type=tipo, default=None)
    a = ap.parse_args(argv)
    if a.selftest:
        from zonas_selftest import rodar
        sys.exit(rodar())
    if not a.entradas:
        ap.error("informe o(s) arquivo(s) do talhão")
    p = Params(**{k: getattr(a, k) for k in ("densidade_ha", "area_min_ha", "area_max_ha", "rigidez_area", "res_m",
                                              "obs_por_janela", "inicio", "fim", "cache_dir") if getattr(a, k) is not None})
    tal = ler_talhoes(zp.listar_arquivos(a.entradas) or a.entradas, a.unir_por_arquivo)
    fonte = FonteSintetica() if a.demo else FonteSTAC()
    erros = []
    res = analisar_fazenda(tal, fonte, p, lambda f, m: print(f"[{f:4.0%}] {m}", flush=True), erros)
    for nome, e in erros:
        print(f"FALHOU: {nome}: {e}")
    if not res:
        sys.exit("Nenhum talhão analisado.")
    Path(a.kml).write_text(kml_zonas(res), encoding="utf-8")
    print(f"KML para upload: {a.kml}")
    if a.camadas:
        Path(a.camadas).write_text(kml_camadas(res), encoding="utf-8")
    if a.painel:
        Path(a.painel).write_bytes(painel_png(res[0]))
    if a.csv:
        tabela_zonas(res).to_csv(a.csv, index=False)
    for r in res:
        s = r.resumo(); t = r.trafego
        print(f"{r.nome}: {s['n_zonas']} zonas, {s['ha_por_amostra']:.1f} ha/amostra, "
              f"{s['fora_dos_limites']} fora dos limites | cenas: {t.ndvi_rede} baixadas, {t.ndvi_cache} do cache, "
              f"{t.cenas_descartadas_nuvem} descartadas por nuvem")


if __name__ == "__main__":
    main()
