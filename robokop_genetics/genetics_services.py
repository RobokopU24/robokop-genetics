from robokop_genetics.services.myvariant import MyVariantService
from robokop_genetics.services.ensembl import EnsemblService
from robokop_genetics.services.hgnc import HGNCService
from robokop_genetics.simple_graph_components import SimpleNode
from robokop_genetics.util import LoggingUtil
from robokop_genetics.genetics_cache import GeneticsCache
from collections import defaultdict
import logging
import os

MYVARIANT = "MyVariant"
ENSEMBL = "Ensembl"

ALL_VARIANT_TO_GENE_SERVICES = [MYVARIANT, ENSEMBL]
BATCHABLE_VARIANT_TO_GENE_SERVES = [MYVARIANT]


class GeneticsServices(object):

    def __init__(self, provided_cache: GeneticsCache = None, use_cache: bool=True, log_file_path: str = None):
        self.logger = LoggingUtil.init_logging(__name__,
                                               logging.INFO,
                                               logFilePath=log_file_path)
        if use_cache:
            if provided_cache:
                self.cache = provided_cache
            else:
                try:
                    self.cache = GeneticsCache(redis_host=os.environ['ROBO_GENETICS_CACHE_HOST'],
                                               redis_port=os.environ['ROBO_GENETICS_CACHE_PORT'],
                                               redis_db=os.environ['ROBO_GENETICS_CACHE_DB'],
                                               redis_password=os.environ['ROBO_GENETICS_CACHE_PASSWORD'],
                                               log_file_path=log_file_path)
                except KeyError:
                    self.logger.debug('ROBO GENETICS CACHE environment variables not set up. No cache activated.')
                    self.cache = None
        else:
            self.cache = None
        self.hgnc = HGNCService(log_file_path)
        self.myvariant = MyVariantService(log_file_path, hgnc_service=self.hgnc)
        self.ensembl = EnsemblService(log_file_path)

    def get_variant_to_gene(self, services: list, variant_nodes: list):
        all_results = defaultdict(list)
        for service in services:
            nodes_that_need_results = []
            if self.cache:
                cache_key = f'{service}_sequence_variant_to_gene'
                cached_results = self.cache.get_service_results(cache_key, [node.id for node in variant_nodes])

                for i, node in enumerate(variant_nodes):
                    cached_result = cached_results[i]
                    if cached_result is not None:
                        all_results[node.id] = cached_result
                    else:
                        nodes_that_need_results.append(node)
            else:
                nodes_that_need_results = variant_nodes

            if service == MYVARIANT:
                variant_dict = {}
                for node in nodes_that_need_results:
                    # TODO this should just pass all the synonyms
                    # but for now the legacy apps with labeled IDs as synonyms wouldn't work
                    # just grabbing plain CURIES that are relevant
                    variant_dict[node.id] = node.get_synonyms_by_prefix('MYVARIANT_HG38')
                new_myvariant_results = self.batch_query_variant_to_gene(MYVARIANT, variant_dict)
                for node_id, results in new_myvariant_results.items():
                    all_results[node_id].extend(results)
                if self.cache:
                    self.cache.set_service_results(MYVARIANT, new_myvariant_results)
            elif service == ENSEMBL:
                new_ensembl_results = {}
                for node in nodes_that_need_results:
                    variant_id = node.id
                    variant_syns = node.get_synonyms_by_prefix('ROBO_VARIANT')
                    new_ensembl_results[variant_id] = self.ensembl.sequence_variant_to_gene(variant_id, variant_syns)
                    all_results[variant_id].extend(new_ensembl_results[variant_id])
                if self.cache:
                    self.cache.set_service_results(ENSEMBL, new_ensembl_results)
        return all_results


    # service: the service to query (from ALL_VARIANT_TO_GENE_SERVICES)
    # variant_id: plain curie string
    # variant_synonyms: a set of synonym curies
    #
    # specify the service and provide variant information to find gene relationships
    # results will be in a list of tuples
    # (edge: SimpleEdge, gene_node: SimpleNode)
    def query_variant_to_gene(self, service: str, variant_id: str, variant_synonyms: set):
        if service == MYVARIANT:
            return self.myvariant.sequence_variant_to_gene(variant_id, variant_synonyms)
        elif service == ENSEMBL:
            return self.ensembl.sequence_variant_to_gene(variant_id, variant_synonyms)
        else:
            self.logger.warning(f'Service ({service}) not found! Variant to gene failed.')

    # variant_dict: a dictionary of variant_id (curie) to variant_synonyms (set of curies)
    # these are the same parameters for get_variant_to_gene
    # returns a dictionary with the variant id curie as keys and the results from get_variant_to_gene as values
    def batch_query_variant_to_gene(self, service: str, variant_dict: dict):
        if service == MYVARIANT:
            return self.myvariant.batch_sequence_variant_to_gene(variant_dict)
        else:
            self.logger.warning(f'Service ({service}) not batch-able! Variant to gene failed.')

    # given a plain string gene_symbol return a valid curie gene ID
    # eg. BRCA1 -> HGNC:1100
    def get_gene_id_from_symbol(self, gene_symbol: str):
        return self.hgnc.get_gene_id_from_symbol(gene_symbol)
