import pytest

from robokop_genetics.genetics_services import *


@pytest.fixture()
def genetics_services():
    return GeneticsServices(use_cache=False)


def test_gene_symbol_to_id(genetics_services):
    gene_id = genetics_services.get_gene_id_from_symbol('ASS1')
    assert gene_id == 'HGNC:758'

    gene_id = genetics_services.get_gene_id_from_symbol('DMD')
    assert gene_id == 'HGNC:2928'

    gene_id = genetics_services.get_gene_id_from_symbol('BRCA1')
    assert gene_id == 'HGNC:1100'

    gene_id = genetics_services.get_gene_id_from_symbol('THISISAFAKEGENE')
    assert gene_id is None
