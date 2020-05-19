import pytest

from robokop_genetics.genetics_normalization import GeneticsNormalizer


@pytest.fixture()
def genetics_normalizer():
    return GeneticsNormalizer()


def test_normalization(genetics_normalizer):
    """Check variant normalization, which mostly occurs through the ClinGen Allele Registry (CAID)"""

    synonym_set = {"CAID:CA128085"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)

    assert normalized_id == 'CAID:CA128085'
    assert normalized_name == 'rs671'
    assert 'HGVS:NC_000012.12:g.111803962G>A' in synonyms
    assert 'CLINVARVARIANT:18390' in synonyms
    assert 'DBSNP:rs671' in synonyms

    synonym_set = {"CLINVARVARIANT:18390"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)

    assert normalized_id == 'CAID:CA128085'
    assert normalized_name == 'rs671'

    # should double check this is what we want - tri-allelic
    synonym_set = {"DBSNP:rs369602258"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)

    # this could be CA6146346 or CA321211
    #assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'MYVARIANT_HG38:chr11:g.68032291C>T' in synonyms
    assert 'MYVARIANT_HG38:chr11:g.68032291C>G' in synonyms
    assert 'CAID:CA6146346' in synonyms
    assert 'CAID:CA321211' in synonyms

    synonym_set = {"HGVS:NC_000023.11:g.32389644G>A"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)
    assert normalized_id == 'CAID:CA267021'
    assert normalized_name == 'rs398123953'
    assert 'MYVARIANT_HG38:chrX:g.32389644G>A' in synonyms
    assert 'CLINVARVARIANT:94623' in synonyms
    assert 'DBSNP:rs398123953' in synonyms
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|A' in synonyms

    synonym_set = {"MYVARIANT_HG19:chr11:g.67799758C>G"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)
    assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'CAID:CA6146346' in synonyms
    assert 'DBSNP:rs369602258' in synonyms
    assert 'HGVS:NC_000011.10:g.68032291C>G' in synonyms
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|G' in synonyms

    synonym_set = {"MYVARIANT_HG38:chr11:g.68032291C>G"}
    normalized_id, normalized_name, synonyms = genetics_normalizer.get_sequence_variant_normalization(synonym_set)
    assert normalized_id == 'CAID:CA6146346'
    assert normalized_name == 'rs369602258'
    assert 'CAID:CA6146346' in synonyms
    assert 'DBSNP:rs369602258' in synonyms
    assert 'HGVS:NC_000011.10:g.68032291C>G' in synonyms
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|G' in synonyms


def test_batch_normalization(genetics_normalizer):

    hgvs_ids = ['NC_000011.10:g.68032291C>G', 'NC_000023.9:g.32317682G>A', 'NC_000017.10:g.43009069G>C',
                'NC_000017.10:g.43009127delG']

    batch_normalizations = genetics_normalizer.get_batch_sequence_variant_normalization(hgvs_ids)

    normalized_id, normalized_name, synonyms = batch_normalizations['HGVS:NC_000023.9:g.32317682G>A']
    assert 'CAID:CA267021' in synonyms
    assert normalized_name == 'rs398123953'

    normalized_id, normalized_name, synonyms = batch_normalizations['HGVS:NC_000011.10:g.68032291C>G']
    assert 'DBSNP:rs369602258' in synonyms
    assert normalized_name == 'rs369602258'

    normalized_id, normalized_name, synonyms = batch_normalizations['HGVS:NC_000017.10:g.43009127delG']
    assert 'DBSNP:rs775219016' in synonyms
    assert 'CAID:CA8609461' in synonyms
    assert 'MYVARIANT_HG38:chr17:g.44931759del' in synonyms
    assert normalized_name == 'rs775219016'
