import pytest
import os
from robokop_genetics.genetics_cache import GeneticsCache
from robokop_genetics.genetics_services import *


@pytest.fixture()
def genetics_cache():
    testing_prefix = 'robokop-genetics-testing-key-'
    if os.environ['ROBO_GENETICS_CACHE_HOST']:
        g_cache = GeneticsCache(
            redis_db=os.environ['ROBO_GENETICS_CACHE_DB'],
            redis_host=os.environ['ROBO_GENETICS_CACHE_HOST'],
            redis_port=os.environ['ROBO_GENETICS_CACHE_PORT'],
            redis_password=os.environ['ROBO_GENETICS_CACHE_PASSWORD'],
            prefix=testing_prefix)
        g_cache.delete_all_keys_with_prefix(testing_prefix)
        return g_cache
    else:
        print('oops, cache environment variables not set!')


@pytest.fixture()
def genetics_services():
    # can use this for troubleshooting
    #log_file_path = './'
    log_file_path = None
    return GeneticsServices(log_file_path)


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

    for node_id in mock_normalizations:
        cached_normalization = genetics_cache.get_normalization(node_id)
        cached_id, cached_name, cached_synonyms = cached_normalization
        expected_id, expected_name, expected_synonyms = mock_normalizations[node_id]
        assert cached_id == expected_id
        assert cached_name == expected_name
        assert cached_synonyms == expected_synonyms


def test_service_results_cache(genetics_cache, genetics_services):

    node_id = 'CAID:CA279509'
    service_key = f'{ENSEMBL}_variant_to_gene'
    robokop_variant_id = f'ROBO_VARIANT:HG38|17|58206171|58206172|A'
    service_results = genetics_services.query_variant_to_gene(ENSEMBL, node_id, {node_id, robokop_variant_id})
    results_dict = {node_id: service_results}
    genetics_cache.set_service_results(service_key, results_dict)

    results_from_cache = genetics_cache.get_service_results(service_key, node_ids=[node_id])
    results = results_from_cache[0]
    for edge, node in results:
        if node.id == "ENSEMBL:ENSG00000108384":
            assert node.name == "RAD51C"
            assert edge.properties['distance'] == 486402
            assert edge.predicate_label == 'nearby_variant_of'

        if node.id == "ENSEMBL:ENSG00000121101":
            assert node.name == "TEX14"
            assert edge.properties['distance'] > 0

    identifiers = [node.id for edge, node in results]
    assert 'ENSEMBL:ENSG00000011143' in identifiers
    assert 'ENSEMBL:ENSG00000121053' in identifiers
    assert 'ENSEMBL:ENSG00000167419' in identifiers
    assert len(identifiers) > 20

