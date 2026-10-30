import io
import zipfile

import pytest

from arquivos_seguros import extrair_zip_seguro


def zip_com(nome, conteudo=b'ok'):
    memoria = io.BytesIO()
    with zipfile.ZipFile(memoria, 'w') as z:
        z.writestr(nome, conteudo)
    memoria.seek(0)
    return memoria


def test_zip_valido(tmp_path):
    arquivos = extrair_zip_seguro(zip_com('pasta/limite.kml'), tmp_path)
    assert len(arquivos) == 1
    assert arquivos[0].read_bytes() == b'ok'


@pytest.mark.parametrize('nome', ['../invasao.txt', '/tmp/invasao.txt', 'C:/invasao.txt', 'pasta/../../invasao.txt'])
def test_bloqueia_caminho_inseguro(tmp_path, nome):
    with pytest.raises(ValueError):
        extrair_zip_seguro(zip_com(nome), tmp_path)


def test_bloqueia_arquivo_muito_grande(tmp_path, monkeypatch):
    import arquivos_seguros
    monkeypatch.setattr(arquivos_seguros, 'MAX_ARQUIVO', 1)
    with pytest.raises(ValueError):
        extrair_zip_seguro(zip_com('teste.txt', b'xx'), tmp_path)
