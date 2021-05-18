import pytest

from robokop_genetics.util import Text
from robokop_genetics.genetics_normalization import GeneticsNormalizer
from robokop_genetics.services.clingen import ClinGenService
import robokop_genetics.node_types as node_types

import requests

"""Check variant synonymization through the ClinGen Allele Registry (CAID)
"""


@pytest.fixture()
def genetics_normalizer():
    return GeneticsNormalizer(use_cache=False)


@pytest.fixture()
def clingen_service():
    return ClinGenService()


def test_errors(clingen_service):

    with pytest.raises(requests.exceptions.RequestException) as re:
        clingen_service.query_service(f'{clingen_service.url}alleles?thisisabrokenrequest!')

    unsupported_ids = ['DBSNP:CA128085', 'DBSNP:rs10791957']
    with pytest.raises(NotImplementedError):
        clingen_service.get_batch_of_synonyms(unsupported_ids)


def test_one_at_a_time_normalization(genetics_normalizer):

    # Some curie types should never be normalized one at a time, these should all fail
    node_id = "CAID:CA128085"
    assert genetics_normalizer.get_sequence_variant_normalization(node_id) == []
    node_id = "HGVS:NC_000023.11:g.32389644G>A"
    assert genetics_normalizer.get_sequence_variant_normalization(node_id) == []
    node_id = "MYVARIANT_HG38:chr11:g.68032291C>G"
    assert genetics_normalizer.get_sequence_variant_normalization(node_id) == []

    node_id = "CLINVARVARIANT:18390"
    normalization_info = genetics_normalizer.get_sequence_variant_normalization(node_id).pop()
    assert normalization_info["id"] == 'CAID:CA128085'
    assert normalization_info["name"] == 'rs671'
    assert 'DBSNP:rs671' in normalization_info["equivalent_identifiers"]

    # rs369602258 is tri-allelic - the following tests show how a specific allele can be normalized from a DBSNP
    # if no allele specified return all CAID and their synonym sets
    node_id = "DBSNP:rs369602258"
    normalizations = genetics_normalizer.get_sequence_variant_normalization(node_id)
    normalized_ids = [norm["id"] for norm in normalizations]
    normalized_names = [norm["name"] for norm in normalizations]
    assert 'CAID:CA6146346' in normalized_ids
    assert 'CAID:CA321211' in normalized_ids
    assert 'rs369602258' in normalized_names

    # if the reference is specified - return all CAID and their synonym sets
    node_id = "DBSNP:rs369602258-C"
    normalizations_2 = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert len(normalizations) == len(normalizations_2)
    normalized_ids = [norm["id"] for norm in normalizations]
    assert 'CAID:CA6146346' in normalized_ids
    assert 'CAID:CA321211' in normalized_ids

    # if a non-reference allele is specified - return only the CAID that matches it
    node_id = "DBSNP:rs369602258-T"
    normalizations = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert len(normalizations) == 1
    assert normalizations[0]["id"] == "CAID:CA321211"
    robo_id = Text.get_curies_by_prefix('ROBO_VARIANT', normalizations[0]["equivalent_identifiers"]).pop()
    assert robo_id.split('|')[-1] == 'T'
    assert robo_id.split('|')[-2] == 'C'

    node_id = "DBSNP:rs369602258-G"
    normalizations = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert len(normalizations) == 1
    assert normalizations[0]["id"] == "CAID:CA6146346"
    robo_id = Text.get_curies_by_prefix('ROBO_VARIANT', normalizations[0]["equivalent_identifiers"]).pop()
    assert robo_id.split('|')[-1] == 'G'
    assert robo_id.split('|')[-2] == 'C'


def test_batch_synonymization(clingen_service):

    hgvs_ids = ['HGVS:NC_000011.10:g.68032291C>G',
                'HGVS:NC_000023.9:g.32317682G>A',
                'HGVS:NC_000017.10:g.43009069G>C',
                'HGVS:NC_000017.10:g.43009127delG',
                'HGVS:NC_000001.40:fakehgvs.1231234A>C']

    batch_synonymizations = clingen_service.get_batch_of_synonyms(hgvs_ids)

    synonyms = batch_synonymizations[0]
    assert 'CAID:CA6146346' in synonyms
    assert 'DBSNP:rs369602258' in synonyms
    assert isinstance(synonyms, set)

    synonyms = batch_synonymizations[1]
    assert 'CAID:CA267021' in synonyms
    assert 'DBSNP:rs398123953' in synonyms
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|G|A' in synonyms

    synonyms = batch_synonymizations[3]
    assert 'CAID:CA8609461' in synonyms
    assert 'DBSNP:rs775219016' in synonyms
    assert 'MYVARIANT_HG38:chr17:g.44931759del' in synonyms

    bad_result = batch_synonymizations[4]
    assert bad_result == set()


def test_batch_normalization(genetics_normalizer):

    hgvs_ids = ['HGVS:NC_000011.10:g.68032291C>G',
                'HGVS:NC_000023.9:g.32317682G>A',
                'HGVS:NC_000017.10:g.43009069G>C',
                'HGVS:NC_000017.10:g.43009127delG',
                'HGVS:NC_000001.40:fakehgvs.1231234A>C']

    batch_normalizations = genetics_normalizer.normalize_variants(hgvs_ids)

    normalization_info = batch_normalizations['HGVS:NC_000011.10:g.68032291C>G'].pop()
    assert normalization_info["id"] == 'CAID:CA6146346'
    assert normalization_info["name"] == 'rs369602258'
    assert node_types.SEQUENCE_VARIANT in normalization_info["type"]
    assert node_types.NAMED_THING in normalization_info["type"]
    assert node_types.MOLECULAR_ENTITY in normalization_info["type"]

    normalization_info = batch_normalizations['HGVS:NC_000023.9:g.32317682G>A'].pop()
    assert normalization_info["id"] == 'CAID:CA267021'
    assert normalization_info["name"] == 'rs398123953'

    normalization_info = batch_normalizations['HGVS:NC_000017.10:g.43009127delG'].pop()
    assert normalization_info["id"] == 'CAID:CA8609461'
    assert normalization_info["name"] == 'rs775219016'
    assert 'MYVARIANT_HG38:chr17:g.44931759del' in normalization_info["equivalent_identifiers"]

    bad_result = batch_normalizations['HGVS:NC_000001.40:fakehgvs.1231234A>C']
    assert bad_result == []


def test_mixed_normalization(genetics_normalizer):

    variant_ids = ['CAID:CA128085',
                   'HGVS:NC_000023.11:g.32389644G>A',
                   'HGVS:NC_000011.10:g.68032291C>T',
                   'HGVS:NC_000011.10:g.68032291C>G',
                   'CLINVARVARIANT:18390',
                   'DBSNP:rs10791957',
                   'BOGUS:rs999999999999']

    normalization_map = genetics_normalizer.normalize_variants(variant_ids)

    assert normalization_map['CAID:CA128085'][0]["id"] == 'CAID:CA128085'
    assert normalization_map['CAID:CA128085'][0]["name"] == 'rs671'
    normalized_synonyms = normalization_map['CAID:CA128085'][0]["equivalent_identifiers"]
    assert 'HGVS:NC_000012.12:g.111803962G>A' in normalized_synonyms
    assert 'CLINVARVARIANT:18390' in normalized_synonyms
    assert 'DBSNP:rs671' in normalized_synonyms
    assert 'MYVARIANT_HG38:chr12:g.111803962G>A' in normalized_synonyms

    assert normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["id"] == 'CAID:CA267021'
    assert normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["name"] == 'rs398123953'
    normalized_synonyms = normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["equivalent_identifiers"]
    assert 'MYVARIANT_HG38:chrX:g.32389644G>A' in normalized_synonyms
    assert 'CLINVARVARIANT:94623' in normalized_synonyms
    assert 'DBSNP:rs398123953' in normalized_synonyms
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|G|A' in normalized_synonyms

    assert normalization_map['HGVS:NC_000011.10:g.68032291C>T'][0]["id"] == "CAID:CA321211"
    assert normalization_map['HGVS:NC_000011.10:g.68032291C>T'][0]["name"] == 'rs369602258'

    assert normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["id"] == 'CAID:CA6146346'
    assert normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["name"] == 'rs369602258'
    normalized_synonyms = normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["equivalent_identifiers"]
    assert 'HGVS:NC_000011.10:g.68032291C>G' in normalized_synonyms
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|C|G' in normalized_synonyms

    assert normalization_map['DBSNP:rs10791957'][0]["id"] == 'CAID:CA15722020'
    normalized_node_types = normalization_map['DBSNP:rs10791957'][0]["type"]
    assert node_types.SEQUENCE_VARIANT in normalized_node_types
    assert node_types.NAMED_THING in normalized_node_types
    assert node_types.MOLECULAR_ENTITY in normalized_node_types
