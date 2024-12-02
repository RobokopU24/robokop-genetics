import pytest

from robokop_genetics.genetics_normalization import GeneticsNormalizer
from robokop_genetics.services.clingen import ClinGenService, ClinGenSynonymizationResult, ClinGenQueryResponse
import robokop_genetics.node_types as node_types


"""Check variant synonymization through the ClinGen Allele Registry (CAID)
"""

@pytest.fixture()
def genetics_normalizer():
    return GeneticsNormalizer(use_cache=False)


@pytest.fixture()
def clingen_service():
    return ClinGenService()


def test_errors(clingen_service):

    clingen_response: ClinGenQueryResponse = clingen_service.query_service(f'{clingen_service.url}alleles?thisisabrokenrequest!')
    assert clingen_response.success is False
    assert clingen_response.error_type == "IncorrectRequest"

    unsupported_ids = ['DBSNP:CA128085', 'DBSNP:rs10791957']
    with pytest.raises(NotImplementedError):
        clingen_service.get_batch_of_synonyms(unsupported_ids)


def test_one_at_a_time_normalization(genetics_normalizer):

    # Some curie types should never be normalized one at a time, these should all fail
    node_id = "CAID:CA128085"
    normalization_result = genetics_normalizer.get_sequence_variant_normalization(node_id).pop()
    assert normalization_result['error_type'] == 'InefficientUsage'
    node_id = "HGVS:NC_000023.11:g.32389644G>A"
    normalization_result = genetics_normalizer.get_sequence_variant_normalization(node_id).pop()
    assert normalization_result['error_type'] == 'InefficientUsage'

    node_id = "CLINVARVARIANT:18390"
    synonymization_results = genetics_normalizer.get_sequence_variant_normalization(node_id)
    found_result = False
    for normalization_info in synonymization_results:
        if 'id' in normalization_info:
            assert normalization_info["id"] == 'CAID:CA128085'
            assert normalization_info["name"] == 'rs671'
            assert 'DBSNP:rs671' in normalization_info["equivalent_identifiers"]
            found_result = True
    assert found_result

    # rs369602258 is tri-allelic - the following tests show how a specific allele can be normalized from a DBSNP
    # if no allele specified return all CAID and their synonym sets
    node_id = "DBSNP:rs369602258"
    normalizations = genetics_normalizer.get_sequence_variant_normalization(node_id)
    normalized_ids = [norm["id"] for norm in normalizations]
    normalized_names = [norm["name"] for norm in normalizations]
    assert 'CAID:CA6146346' in normalized_ids
    assert 'CAID:CA321211' in normalized_ids
    assert 'rs369602258' in normalized_names

    # if the allele preference is not found as an alt allele (in this case it's the reference)
    # return all CAID and their synonym sets
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

    robo_id = normalizations[0]["robokop_variant_id"]
    assert robo_id.split('|')[-1] == 'T'
    assert robo_id.split('|')[-2] == 'C'

    node_id = "DBSNP:rs369602258-G"
    normalizations = genetics_normalizer.get_sequence_variant_normalization(node_id)
    assert len(normalizations) == 1
    assert normalizations[0]["id"] == "CAID:CA6146346"
    robo_id = normalizations[0]["robokop_variant_id"]
    assert robo_id.split('|')[-1] == 'G'
    assert robo_id.split('|')[-2] == 'C'


def test_batch_synonymization(clingen_service):

    hgvs_ids = ['HGVS:NC_000011.10:g.68032291C>G',
                'HGVS:NC_000023.9:g.32317682G>A',
                'HGVS:NC_000017.10:g.43009069G>C',
                'HGVS:NC_000017.10:g.43009127delG',
                'HGVS:NC_000001.40:fakehgvs.1231234A>C']

    batch_synonymizations = clingen_service.get_batch_of_synonyms(hgvs_ids)

    synonymization_result: ClinGenSynonymizationResult = batch_synonymizations[0]
    assert synonymization_result.success
    assert 'CAID:CA6146346' == synonymization_result.id
    assert 'CAID:CA6146346' not in synonymization_result.equivalent_identifiers
    assert 'DBSNP:rs369602258' in synonymization_result.equivalent_identifiers
    assert isinstance(synonymization_result.equivalent_identifiers, list)

    synonymization_result: ClinGenSynonymizationResult = batch_synonymizations[1]
    assert 'CAID:CA267021' == synonymization_result.id
    assert 'DBSNP:rs398123953' in synonymization_result.equivalent_identifiers
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|G|A' == synonymization_result.robokop_variant_id

    synonymization_result: ClinGenSynonymizationResult = batch_synonymizations[3]
    assert 'CAID:CA8609461' == synonymization_result.id
    assert 'DBSNP:rs775219016' in synonymization_result.equivalent_identifiers

    synonymization_result: ClinGenSynonymizationResult = batch_synonymizations[4]
    assert synonymization_result.success is False
    assert synonymization_result.error_type == 'HgvsParsingError'


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
    assert node_types.SEQUENCE_VARIANT in normalization_info["category"]
    assert node_types.NAMED_THING in normalization_info["category"]
    assert node_types.BIOLOGICAL_ENTITY in normalization_info["category"]

    normalization_info = batch_normalizations['HGVS:NC_000023.9:g.32317682G>A'].pop()
    assert normalization_info["id"] == 'CAID:CA267021'
    assert normalization_info["name"] == 'rs398123953'

    normalization_info = batch_normalizations['HGVS:NC_000017.10:g.43009127delG'].pop()
    assert normalization_info["id"] == 'CAID:CA8609461'
    assert normalization_info["name"] == 'rs775219016'

    normalization_info = batch_normalizations['HGVS:NC_000001.40:fakehgvs.1231234A>C'].pop()
    assert 'error_type' in normalization_info
    assert normalization_info['error_type'] == 'HgvsParsingError'


def test_mixed_normalization(genetics_normalizer):

    variant_ids = ['CAID:CA128085',
                   'HGVS:NC_000023.11:g.32389644G>A',
                   'HGVS:NC_000011.10:g.68032291C>T',
                   'HGVS:NC_000011.10:g.68032291C>G',
                   'CLINVARVARIANT:18390',
                   'DBSNP:rs10791957',
                   'BOGUS:rs999999999999',
                   'DBSNP:rs199745043-AG',
                   'DBSNP:rs3180018']

    normalization_map = genetics_normalizer.normalize_variants(variant_ids)

    assert normalization_map['CAID:CA128085'][0]["id"] == 'CAID:CA128085'
    assert normalization_map['CAID:CA128085'][0]["name"] == 'rs671'
    print(normalization_map)
    equivalent_identifiers = normalization_map['CAID:CA128085'][0]["equivalent_identifiers"]
    assert 'CLINVARVARIANT:18390' in equivalent_identifiers
    assert 'DBSNP:rs671' in equivalent_identifiers
    hgvs_ids = normalization_map['CAID:CA128085'][0]["hgvs"]
    assert 'HGVS:NC_000012.12:g.111803962G>A' in hgvs_ids

    assert normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["id"] == 'CAID:CA267021'
    assert normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["name"] == 'rs398123953'
    equivalent_identifiers = normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["equivalent_identifiers"]
    assert 'CLINVARVARIANT:94623' in equivalent_identifiers
    assert 'DBSNP:rs398123953' in equivalent_identifiers
    robokop_variant_id = normalization_map['HGVS:NC_000023.11:g.32389644G>A'][0]["robokop_variant_id"]
    assert 'ROBO_VARIANT:HG38|X|32389643|32389644|G|A' == robokop_variant_id

    assert normalization_map['HGVS:NC_000011.10:g.68032291C>T'][0]["id"] == "CAID:CA321211"
    assert normalization_map['HGVS:NC_000011.10:g.68032291C>T'][0]["name"] == 'rs369602258'

    assert normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["id"] == 'CAID:CA6146346'
    assert normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["name"] == 'rs369602258'
    hgvs_ids = normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["hgvs"]
    assert 'HGVS:NC_000011.10:g.68032291C>G' in hgvs_ids
    robokop_variant_id = normalization_map['HGVS:NC_000011.10:g.68032291C>G'][0]["robokop_variant_id"]
    assert 'ROBO_VARIANT:HG38|11|68032290|68032291|C|G' == robokop_variant_id

    assert normalization_map['DBSNP:rs10791957'][0]["id"] == 'CAID:CA1980501278'
    normalized_node_types = normalization_map['DBSNP:rs10791957'][0]["category"]
    assert node_types.SEQUENCE_VARIANT in normalized_node_types
    assert node_types.NAMED_THING in normalized_node_types
    assert node_types.BIOLOGICAL_ENTITY in normalized_node_types
    assert normalization_map['DBSNP:rs10791957'][1]["id"] == 'CAID:CA15722020'

    assert normalization_map['DBSNP:rs199745043-AG'][0]["id"] == 'CAID:CA101278074'
    assert len(normalization_map['DBSNP:rs199745043-AG']) == 1

    assert normalization_map['BOGUS:rs999999999999'][0]["error_type"] == 'UnsupportedPrefix'
    assert 'BOGUS' in normalization_map['BOGUS:rs999999999999'][0]["error_message"]

    assert normalization_map['DBSNP:rs3180018'][0]["error_type"] == 'NotFound'
