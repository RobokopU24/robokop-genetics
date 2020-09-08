import pytest
import os
from robokop_genetics.genetics_cache import GeneticsCache
from robokop_genetics.genetics_services import *


@pytest.fixture()
def genetics_cache():
    testing_prefix = 'robo-testing-key-'
    if 'ROBO_GENETICS_CACHE_HOST' in os.environ:
        testing_cache = GeneticsCache(prefix=testing_prefix)
        testing_cache.delete_all_keys_with_prefix(testing_prefix)
        return testing_cache
    else:
        pytest.fail('Cache environment variables not set! Cache can not be utilized.')


@pytest.fixture()
def genetics_services():
    return GeneticsServices()


@pytest.fixture()
def mock_normalizations():
    mock_normalizations = dict()
    mock_normalizations['TESTINGCURIE:10'] = ('TESTINGCURIE:10',
                                             'TESTING NAME 1',
                                             ['TESTINGCURIE:10',
                                              'TESTINGCURIE:11',
                                              'TESTINGCURIE:12'])
    mock_normalizations['TESTINGCURIE:20'] = ('TESTINGCURIE:21',
                                             'TESTING NAME 2',
                                             ['TESTINGCURIE:20',
                                              'TESTINGCURIE:21',
                                              'TESTINGCURIE:22'])
    mock_normalizations['TESTINGCURIE:30'] = ('TESTINGCURIE:33',
                                             'TESTING NAME 3',
                                             ['TESTINGCURIE:30',
                                              'TESTINGCURIE:31',
                                              'TESTINGCURIE:33'])
    mock_normalizations['TESTINGCURIE:40'] = ('TESTINGCURIE:45',
                                             'TESTING NAME 4',
                                             ['TESTINGCURIE:40',
                                              'TESTINGCURIE:41',
                                              'TESTINGCURIE:42',
                                              'TESTINGCURIE:45'])
    return mock_normalizations


def test_normalization_cache(genetics_cache, mock_normalizations):
    node_id = next(iter(mock_normalizations.keys()))
    cached_normalization = genetics_cache.get_normalization(node_id)
    assert cached_normalization is None

    genetics_cache.set_normalization(node_id, mock_normalizations[node_id])
    cached_normalization = genetics_cache.get_normalization(node_id)
    cached_id, cached_name, cached_synonyms = cached_normalization
    expected_id, expected_name, expected_synonyms = mock_normalizations[node_id]
    assert cached_id == expected_id
    assert cached_name == expected_name
    assert cached_synonyms == expected_synonyms


def test_batch_normalization_cache(genetics_cache, mock_normalizations):
    for node_id in mock_normalizations:
        cached_normalization = genetics_cache.get_normalization(node_id)
        assert cached_normalization is None

    genetics_cache.set_batch_normalization(mock_normalizations)

    list_of_node_ids = list(mock_normalizations.keys())
    batch_normalizations = genetics_cache.get_batch_normalization(list_of_node_ids)

    for i, node_id in enumerate(list_of_node_ids):
        cached_id, cached_name, cached_synonyms = batch_normalizations[i]
        expected_id, expected_name, expected_synonyms = mock_normalizations[node_id]
        assert cached_id == expected_id
        assert cached_name == expected_name
        assert cached_synonyms == expected_synonyms


def test_service_results_cache(genetics_cache, genetics_services):

    # TODO these service queries should be mocked instead of calling the service
    results_dict = {}
    node_id = 'CAID:CA279509'
    robokop_variant_id = f'ROBO_VARIANT:HG38|17|58206171|58206172|A'
    service_results = genetics_services.query_variant_to_gene(ENSEMBL, node_id, {node_id, robokop_variant_id})
    results_dict[node_id] = service_results

    node_id_2 = 'FAKECURIE:39'
    robokop_variant_id = f'ROBO_VARIANT:HG38|X|32389643|32389644|A'
    service_results = genetics_services.query_variant_to_gene(ENSEMBL, node_id_2, {node_id_2, robokop_variant_id})
    results_dict[node_id_2] = service_results

    service_key = f'{ENSEMBL}_variant_to_gene'
    genetics_cache.set_service_results(service_key, results_dict)

    node_id_3 = 'MADEUP:1000'  # this shouldnt be in the cache
    results_from_cache = genetics_cache.get_service_results(service_key, node_ids=[node_id, node_id_2, node_id_3])
    results = results_from_cache[0]
    found_1 = False
    found_2 = False
    for edge, node in results:
        if node.id == "ENSEMBL:ENSG00000108384":
            assert node.name == "RAD51C"
            assert edge.properties['distance'] == 486402
            assert edge.predicate_label == 'upstream_gene_variant'
            found_1 = True

        if node.id == "ENSEMBL:ENSG00000121101":
            assert node.name == "TEX14"
            assert edge.properties['distance'] > 0
            found_2 = True
    assert found_1 and found_2

    identifiers = [node.id for edge, node in results]
    assert 'ENSEMBL:ENSG00000011143' in identifiers
    assert 'ENSEMBL:ENSG00000121053' in identifiers
    assert 'ENSEMBL:ENSG00000167419' in identifiers
    assert len(identifiers) > 20

    results = results_from_cache[1]
    identifiers = [node.id for edge, node in results]
    assert 'ENSEMBL:ENSG00000198947' in identifiers

    # this one had no robo variant key so it shouldnt have a result
    results = results_from_cache[2]
    assert results is None

    node_id_3 = 'CAID:CA6451230'
    service_key = f'{MYVARIANT}_variant_to_gene'
    myvariant_curie = "MYVARIANT_HG38:chr12:g.11091595T>C"
    service_results = genetics_services.query_variant_to_gene(MYVARIANT, node_id_3, {node_id_3, myvariant_curie})
    identifiers = [node.id for edge, node in service_results]
    assert 'HGNC:18875' in identifiers
    predicates = [edge.predicate_id for edge, node in service_results]
    assert 'SNPEFF:missense_variant' in predicates
    results_dict = {node_id_3: service_results}
    genetics_cache.set_service_results(service_key, results_dict)

    results_from_cache = genetics_cache.get_service_results(service_key, node_ids=[node_id_3])
    results = results_from_cache[0]
    identifiers = [node.id for edge, node in results]
    assert 'HGNC:9366' in identifiers
    predicates = [edge.predicate_id for edge, node in results]
    assert 'SNPEFF:intron_variant' in predicates
