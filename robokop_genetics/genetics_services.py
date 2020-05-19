from robokop_genetics.services.myvariant import MyVariantService
from robokop_genetics.services.ensembl import EnsemblService
from robokop_genetics.services.hgnc import HGNCService
from robokop_genetics.util import LoggingUtil
import logging

MYVARIANT = "MyVariant"
ENSEMBL = "Ensembl"

ALL_VARIANT_TO_GENE_SERVICES = [MYVARIANT, ENSEMBL]


class GeneticsServices(object):

    def __init__(self, log_file_path: str = None):
        if log_file_path is None:
            self.logging_on = False
        else:
            self.logging_on = True
            self.logger = LoggingUtil.init_logging(__name__,
                                                   logging.INFO,
                                                   logFilePath=log_file_path)
        self.hgnc = HGNCService(log_file_path)
        self.myvariant = MyVariantService(log_file_path, hgnc_service=self.hgnc)
        self.ensembl = EnsemblService(log_file_path)

    # service: the service to query (from ALL_VARIANT_TO_GENE_SERVICES)
    # variant_id: plain curie string
    # variant_synonyms: a set of synonym curies
    #
    # specify the service and provide variant information to find gene relationships
    # results will be in a list of tuples
    # (edge: SimpleEdge, gene_node: SimpleNode)
    def get_variant_to_gene(self, service: str, variant_id: str, variant_synonyms: set):
        if service == MYVARIANT:
            return self.myvariant.sequence_variant_to_gene(variant_id, variant_synonyms)
        elif service == ENSEMBL:
            return self.ensembl.sequence_variant_to_gene(variant_id, variant_synonyms)
        else:
            if self.logging_on:
                self.logger.warning(f'Service ({service}) not found! Variant to gene failed.')

    # variant_dict: a dictionary of variant_id (curie) to variant_synonyms (set of curies)
    # these are the same parameters for get_variant_to_gene
    # returns a dictionary with the variant id curie as keys and the results from get_variant_to_gene as values
    def batch_variant_to_gene(self, service: str, variant_dict: dict):
        if service == MYVARIANT:
            return self.myvariant.batch_sequence_variant_to_gene(variant_dict)
        else:
            if self.logging_on:
                self.logger.warning(f'Service ({service}) not batch-able! Variant to gene failed.')

    # given a plain string gene_symbol return a valid curie gene ID
    # eg. BRCA1 -> HGNC:1100
    def get_gene_id_from_symbol(self, gene_symbol: str):
        return self.hgnc.get_gene_id_from_symbol(gene_symbol)
