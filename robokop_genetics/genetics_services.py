from robokop_genetics.services.ensembl import EnsemblService
from robokop_genetics.services.hgnc import HGNCService
from robokop_genetics.util import LoggingUtil
from robokop_genetics.genetics_cache import GeneticsCache
from collections import defaultdict
import logging


ENSEMBL = "Ensembl"


class GeneticsServices(object):

    logger = LoggingUtil.init_logging(__name__,
                                      logging.INFO,
                                      log_file_path=LoggingUtil.get_logging_path())

    def __init__(self, use_cache: bool = True):

        if use_cache:
            self.cache = GeneticsCache()
            self.logger.info('Robokop Genetics Services initialized with cache activated.')
        else:
            self.cache = None
            self.logger.info('Robokop Genetics Services initialized with no cache activated.')

        self.hgnc = HGNCService()
        self.ensembl = EnsemblService(temp_dir=LoggingUtil.get_logging_path())

    def get_variant_to_gene(self, services: list, variant_nodes: list):
        self.logger.info(f'Get variant to gene called on {len(variant_nodes)} nodes.')
        all_results = defaultdict(list)
        for service in services:
            if self.cache:
                cache_key = f'{service}_sequence_variant_to_gene'
                cached_results = self.cache.get_service_results(cache_key, [node.id for node in variant_nodes])

                nodes_that_need_results = []
                for i, node in enumerate(variant_nodes):
                    cached_result = cached_results[i]
                    if cached_result is not None:
                        all_results[node.id].extend(cached_result)
                    else:
                        nodes_that_need_results.append(node)
                self.logger.info(f'{service} variant to gene found results for {len(variant_nodes) - len(nodes_that_need_results)} nodes in the cache.')
            else:
                nodes_that_need_results = variant_nodes

            if service == ENSEMBL:
                new_ensembl_results = {}
                counter = 0
                for node in nodes_that_need_results:
                    variant_id = node.id
                    variant_syns = node.get_synonyms_by_prefix('ROBO_VARIANT')
                    new_ensembl_results[variant_id] = self.ensembl.sequence_variant_to_gene(variant_id, variant_syns)
                    all_results[variant_id].extend(new_ensembl_results[variant_id])
                    counter += 1
                    if counter == 10000 and self.cache:
                        self.cache.set_service_results(cache_key, new_ensembl_results)
                        new_ensembl_results = {}
                        counter = 0

                if counter > 0 and self.cache:
                    self.cache.set_service_results(cache_key, new_ensembl_results)

        return all_results

    # service: the service to query (from ALL_VARIANT_TO_GENE_SERVICES)
    # variant_id: plain curie string
    # variant_synonyms: a set of synonym curies
    #
    # specify the service and provide variant information to find gene relationships
    # results will be in a list of tuples
    # (edge: SimpleEdge, gene_node: SimpleNode)
    def query_variant_to_gene(self, service: str, variant_id: str, variant_synonyms: set):
        if service == ENSEMBL:
            return self.ensembl.sequence_variant_to_gene(variant_id, variant_synonyms)
        else:
            self.logger.warning(f'Service ({service}) not found! Variant to gene failed.')

    # given a plain string gene_symbol return a valid curie gene ID
    # eg. BRCA1 -> HGNC:1100
    def get_gene_id_from_symbol(self, gene_symbol: str):
        return self.hgnc.get_gene_id_from_symbol(gene_symbol)
