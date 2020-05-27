import pytest

from robokop_genetics.genetics_normalization import GeneticsNormalizer
from robokop_genetics.simple_graph_components import SimpleNode
import robokop_genetics.node_types as node_types


@pytest.fixture()
def genetics_normalizer():
    return GeneticsNormalizer(use_cache=False)


def test_normalization(genetics_normalizer):
    """Check variant normalization, which mostly occurs through the ClinGen Allele Registry (CAID)"""

    node_id = "CAID:CA128085"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)

    assert normalized_id == 'CAID:CA128085'
    assert normalized_name == 'rs671'
    assert 'HGVS:NC_000012.12:g.111803962G>A' in synonyms
    assert 'CLINVARVARIANT:18390' in synonyms
    assert 'DBSNP:rs671' in synonyms

    node_id = "CLINVARVARIANT:18390"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)

    assert normalized_id == 'CAID:CA128085'
    assert normalized_name == 'rs671'

    # should double check this is what we want - tri-allelic
    node_id = "DBSNP:rs369602258"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)

    # this could be CA6146346 or CA321211
    #assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'MYVARIANT_HG38:chr11:g.68032291C>T' in synonyms
    assert 'MYVARIANT_HG38:chr11:g.68032291C>G' in synonyms
    assert 'CAID:CA6146346' in synonyms
    assert 'CAID:CA321211' in synonyms

    node_id = "HGVS:NC_000023.11:g.32389644G>A"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert normalized_id == 'CAID:CA267021'
    assert normalized_name == 'rs398123953'
    assert 'MYVARIANT_HG38:chrX:g.32389644G>A' in synonyms
    assert 'CLINVARVARIANT:94623' in synonyms
    assert 'DBSNP:rs398123953' in synonyms
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|A' in synonyms

    node_id = "MYVARIANT_HG19:chr11:g.67799758C>G"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'CAID:CA6146346' in synonyms
    assert 'DBSNP:rs369602258' in synonyms
    assert 'HGVS:NC_000011.10:g.68032291C>G' in synonyms
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|G' in synonyms

    node_id = "MYVARIANT_HG38:chr11:g.68032291C>G"
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'CAID:CA6146346' in synonyms
    assert 'DBSNP:rs369602258' in synonyms
    assert 'HGVS:NC_000011.10:g.68032291C>G' in synonyms
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|G' in synonyms


def test_batch_normalization(genetics_normalizer):

    hgvs_ids = ['NC_000011.10:g.68032291C>G',
                'NC_000023.9:g.32317682G>A',
                'NC_000017.10:g.43009069G>C',
                'NC_000017.10:g.43009127delG',
                'NC_000001.40:fakehgvs.1231234A>C']

    batch_normalizations = genetics_normalizer.get_batch_sequence_variant_normalization(hgvs_ids)

    normalized_id, normalized_name, synonyms = batch_normalizations[0]
    assert 'DBSNP:rs369602258' in synonyms
    assert normalized_name == 'rs369602258'
    assert isinstance(synonyms, list)

    normalized_id, normalized_name, synonyms = batch_normalizations[1]
    assert 'CAID:CA267021' in synonyms
    assert normalized_name == 'rs398123953'

    normalized_id, normalized_name, synonyms = batch_normalizations[3]
    assert 'DBSNP:rs775219016' in synonyms
    assert 'CAID:CA8609461' in synonyms
    assert 'MYVARIANT_HG38:chr17:g.44931759del' in synonyms
    assert normalized_name == 'rs775219016'

    normalized_id, normalized_name, synonyms = batch_normalizations[4]
    assert normalized_id == 'HGVS:NC_000001.40:fakehgvs.1231234A>C'
    assert normalized_name == 'NC_000001.40:fakehgvs.1231234A>C'
    assert len(synonyms) == 1

def test_node_based_normalization(genetics_normalizer):
    node = SimpleNode('CAID:CA128085', node_types.SEQUENCE_VARIANT, 'CA128085')
    genetics_normalizer.normalize(node)

    assert node.id == 'CAID:CA128085'
    assert node.name == 'rs671'
    assert len(node.synonyms) > 4
    assert 'HGVS:NC_000012.12:g.111803962G>A' in node.synonyms

    node = SimpleNode('HGVS:NC_000023.11:g.32389644G>A', node_types.SEQUENCE_VARIANT, 'FakeName')
    genetics_normalizer.normalize(node)

    assert node.id == 'CAID:CA267021'
    assert node.name == 'rs398123953'
    assert len(node.synonyms) > 6


def test_node_based_batch_normalization(genetics_normalizer):
    nodes = list()
    nodes.append(SimpleNode('CAID:CA128085', node_types.SEQUENCE_VARIANT, 'CA128085'))
    nodes.append(SimpleNode('HGVS:NC_000023.11:g.32389644G>A', node_types.SEQUENCE_VARIANT, 'FakeName'))
    nodes.append(SimpleNode('HGVS:NC_000011.10:g.68032291C>G', node_types.SEQUENCE_VARIANT, 'FakeName2'))
    nodes.append(SimpleNode('CLINVARVARIANT:18390', node_types.SEQUENCE_VARIANT, 'FakeName3'))
    nodes.append(SimpleNode('DBSNP:rs10791957', node_types.SEQUENCE_VARIANT, 'FakeName4'))
    nodes.append(SimpleNode('BOGUS:rs999999999999', node_types.SEQUENCE_VARIANT, 'BogusName1'))
    genetics_normalizer.batch_normalize(nodes)

    node = nodes[0]
    assert node.id == 'CAID:CA128085'
    assert node.name == 'rs671'
    assert 'HGVS:NC_000012.12:g.111803962G>A' in node.synonyms

    node = nodes[1]
    assert node.id == 'CAID:CA267021'
    assert node.name == 'rs398123953'

    node = nodes[4]
    assert node.id == 'CAID:CA15722020'
    assert node.name == 'rs10791957'
    assert 'HGVS:NC_000011.10:g.68100081C>A' in node.synonyms

    node = nodes[5]
    assert node.id == 'BOGUS:rs999999999999'
    assert node.name == 'rs999999999999'