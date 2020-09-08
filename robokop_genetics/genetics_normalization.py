from robokop_genetics.services.clingen import ClinGenService
from robokop_genetics.genetics_cache import GeneticsCache
from robokop_genetics.simple_graph_components import SimpleNode
from robokop_genetics.util import LoggingUtil, Text
from typing import Set, List
import logging


class GeneticsNormalizer(object):

    def __init__(self, use_cache: bool=True):
        log_file_path = LoggingUtil.get_logging_path()
        self.logger = LoggingUtil.init_logging(__name__,
                                               logging.INFO,
                                               log_file_path=log_file_path)
        if use_cache:
            self.cache = GeneticsCache()
            self.logger.info('Robokop Genetics Normalizer initialized with cache activated.')
        else:
            self.cache = None
            self.logger.info('Robokop Genetics Normalizer initialized with no cache activated.')

        self.clingen = ClinGenService()

    def normalize(self, node: SimpleNode):
        normalization = self.cache.get_normalization(node.id) if self.cache else None
        if normalization is None:
            normalization = self.get_sequence_variant_normalization(node.id)
            if self.cache:
                self.cache.set_normalization(node.id, normalization)
        self.apply_normalization(node, normalization)

    def batch_normalize(self, nodes: list):
        node_ids = [node.id for node in nodes]
        nodes_for_batch_normalizing = []
        hgvs_list_for_batch_normalizing = []
        new_normalizations = {}
        cached_normalizations = self.cache.get_batch_normalization(node_ids) if self.cache else None
        cached_result_count = 0
        apply_normalization = self.apply_normalization
        for i, current_node in enumerate(nodes):
            normalization = cached_normalizations[i] if cached_normalizations else None
            if normalization is None:
                hgvs_curies = current_node.get_synonyms_by_prefix('HGVS')
                if hgvs_curies:
                    nodes_for_batch_normalizing.append(current_node)
                    hgvs_list_for_batch_normalizing.append(Text.un_curie(next(iter(hgvs_curies))))
                else:
                    # this node has no HGVS id so it cant be batched, go ahead and do it by itself
                    normalization = self.get_sequence_variant_normalization(current_node.id)
                    new_normalizations[current_node.id] = normalization
            else:
                cached_result_count += 1
            if normalization is not None:
                apply_normalization(current_node, normalization)
        self.logger.info(f'Batch normalizing found {cached_result_count}/{len(nodes)} nodes in the cache.')
        batch_normalizations = self.get_batch_sequence_variant_normalization(hgvs_list_for_batch_normalizing)
        for i, normalization in enumerate(batch_normalizations):
            current_node = nodes_for_batch_normalizing[i]
            if normalization is None:
                if len(current_node.synonyms) > 1:
                    # we only get here if the node had an HGVS but didnt normalize using that
                    normalization = self.get_sequence_variant_normalization(current_node.id)
                else:
                    # give up and accept the node as is,
                    # it only has one synonym and hasn't been successfully normalized
                    normalization = current_node.id, Text.un_curie(current_node.id), current_node.synonyms
            new_normalizations[current_node.id] = normalization
            apply_normalization(current_node, normalization)
            # this actually does something, id could have changed in apply_normalization, we want both keys
            new_normalizations[current_node.id] = normalization
        if self.cache:
            self.cache.set_batch_normalization(new_normalizations)

    def apply_normalization(self, node: SimpleNode, normalization: tuple):
        normalized_id, normalized_name, synonyms = normalization
        node.id = normalized_id
        node.name = normalized_name
        node.add_synonyms(set(synonyms))

    # node_id: the id of the node that needs normalizing
    #
    # returns a normalization tuple of:
    # normalized_id - string (a curie)
    # normalized_name - string
    # normalized_synonyms - a set of curie strings
    def get_sequence_variant_normalization(self, node_id: str):
        normalized_synonyms = self.update_sequence_variant_synonyms({node_id})
        normalized_id, normalized_name = self.get_id_and_name_from_synonyms(normalized_synonyms)
        return normalized_id, normalized_name, list(normalized_synonyms)

    # Given a list of plain hgvs ids,
    # return a list of corresponding normalization tuples as values
    # (normalized_id, normalized_name, normalized_synonyms)
    def get_batch_sequence_variant_normalization(self, hgvs_list: List):
        normalization_list = []
        synonym_list = self.clingen.get_batch_of_synonyms(hgvs_list)
        for i, normalized_synonyms in enumerate(synonym_list):
            if normalized_synonyms:
                normalized_id, normalized_name = self.get_id_and_name_from_synonyms(normalized_synonyms)
            else:
                normalized_id = f'HGVS:{hgvs_list[i]}'
                normalized_name = hgvs_list[i]
                normalized_synonyms = {normalized_id}
            normalization_list.append((normalized_id, normalized_name, list(normalized_synonyms)))
        return normalization_list

    # helper function to split clingen calls based on available curie types
    def update_sequence_variant_synonyms(self, synonyms: set):
        caid_curies = Text.get_curies_by_prefix('CAID', synonyms)
        if caid_curies:
            synonyms.update(self.clingen.get_synonyms_by_caid(Text.un_curie(caid_curies.pop())))
        else:
            synonyms.update(self.clingen.get_synonyms_by_other_ids(synonyms))
        return synonyms

    # extract the preferred curie and name from the synonym set
    # we prefer CAID for the ID and DBSNP as the name if available
    def get_id_and_name_from_synonyms(self, synonyms: Set):
        normalized_id = None
        caid_curies = Text.get_curies_by_prefix('CAID', synonyms)
        if caid_curies:
            caid_curie = caid_curies.pop()
            normalized_id = caid_curie
            normalized_name = Text.un_curie(caid_curie)

        rsid_curies = Text.get_curies_by_prefix('DBSNP', synonyms)
        if rsid_curies:
            rsid_curie = rsid_curies.pop()
            normalized_name = Text.un_curie(rsid_curie)
            if not normalized_id:
                normalized_id = rsid_curie

        if not normalized_id:
            # we didn't find a CAID or rsid, just take the first one as an arbitrary id/name
            arbitrary_syn = next(iter(synonyms))
            normalized_id = arbitrary_syn
            normalized_name = Text.un_curie(arbitrary_syn)

        return normalized_id, normalized_name
