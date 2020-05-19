from robokop_genetics.services.clingen import ClinGenService
from robokop_genetics.util import Text
from typing import Set, List


class GeneticsNormalizer(object):

    def __init__(self, log_file_path: str = None):
        self.clingen = ClinGenService(log_file_path)


    # provided_synonyms: a set of synonym curies for a sequence variant
    #
    # returns a tuple of:
    # normalized_id - string (a curie)
    # normalized_name - string
    # normalized_synonyms - a set of curie strings
    def get_sequence_variant_normalization(self, provided_synonyms: Set):
        normalized_synonyms = self.update_sequence_variant_synonyms(provided_synonyms)
        normalized_id, normalized_name = self.get_id_and_name_from_synonyms(normalized_synonyms)
        return normalized_id, normalized_name, normalized_synonyms

    # Given a list of plain hgvs ids,
    # return a dictionary with hgvs curies as keys
    # and corresponding normalization tuples as values
    # (normalized_id, normalized_name, normalized_synonyms)
    def get_batch_sequence_variant_normalization(self, hgvs_list: List):
        normalization_dict = {}
        synonym_dict = self.clingen.get_batch_of_synonyms(hgvs_list)
        for curie, normalized_synonyms in synonym_dict.items():
            normalized_id, normalized_name = self.get_id_and_name_from_synonyms(normalized_synonyms)
            normalization_dict[curie] = normalized_id, normalized_name, normalized_synonyms
        return normalization_dict

    # helper function to split clingen calls based on available curie types
    def update_sequence_variant_synonyms(self, synonyms: Set):
        caid_curies = Text.get_curies_by_prefix('CAID', synonyms)
        if caid_curies:
            synonyms.update(self.clingen.get_synonyms_by_caid(Text.un_curie(caid_curies.pop())))
        else:
            synonyms.update(self.clingen.get_synonyms_by_other_ids(synonyms))
        return synonyms

    # extract the preferred curie and name from the synonym set
    # we prefer CAID for the ID and DBSNP as the name if available
    def get_id_and_name_from_synonyms(self, synonyms: Set):
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
