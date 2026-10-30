"""Extração controlada de ZIP/KMZ enviados pelo usuário.

Não extrai links, caminhos absolutos, travessia de diretórios ou arquivos
excessivos. Os limites protegem especialmente aplicativos públicos.
"""
from __future__ import annotations

import stat
import zipfile
from pathlib import Path, PurePosixPath

MAX_MEMBROS = 500
MAX_ARQUIVO = 150 * 1024 * 1024
MAX_TOTAL = 500 * 1024 * 1024


def extrair_zip_seguro(arquivo: str | Path, destino: str | Path) -> list[Path]:
    destino = Path(destino).resolve()
    destino.mkdir(parents=True, exist_ok=True)
    extraidos: list[Path] = []
    with zipfile.ZipFile(arquivo) as z:
        membros = z.infolist()
        if len(membros) > MAX_MEMBROS:
            raise ValueError(f"Arquivo compactado excede {MAX_MEMBROS} entradas.")
        total = 0
        for membro in membros:
            nome = membro.filename.replace('\\', '/')
            caminho = PurePosixPath(nome)
            if (not nome or nome.startswith('/') or '..' in caminho.parts or
                    ':' in caminho.parts[0] or '\x00' in nome):
                raise ValueError(f"Caminho inseguro no ZIP: {nome!r}")
            modo = (membro.external_attr >> 16)
            if stat.S_ISLNK(modo):
                raise ValueError(f"Links simbólicos não são aceitos: {nome!r}")
            if membro.is_dir():
                continue
            if membro.file_size > MAX_ARQUIVO:
                raise ValueError(f"Arquivo interno excede {MAX_ARQUIVO // (1024*1024)} MB: {nome!r}")
            total += membro.file_size
            if total > MAX_TOTAL:
                raise ValueError("Conteúdo descompactado excede 500 MB.")
            alvo = (destino / nome).resolve()
            if not alvo.is_relative_to(destino):
                raise ValueError(f"Caminho fora do destino: {nome!r}")
            alvo.parent.mkdir(parents=True, exist_ok=True)
            # Cópia por blocos: checa o tamanho real em vez de confiar só no cabeçalho.
            escritos = 0
            with z.open(membro) as origem, alvo.open('wb') as saida:
                while bloco := origem.read(1024 * 1024):
                    escritos += len(bloco)
                    if escritos > MAX_ARQUIVO or escritos > membro.file_size:
                        alvo.unlink(missing_ok=True)
                        raise ValueError(f"Tamanho inconsistente no ZIP: {nome!r}")
                    saida.write(bloco)
            extraidos.append(alvo)
    return extraidos
