import robokop_genetics.node_types as node_types
from robokop_genetics.simple_graph_components import SimpleNode
from robokop_genetics.util import Text, LoggingUtil
from math import ceil
import logging
import requests

# other classes should check this list before calling get_batch_of_synonyms
batchable_variant_curie_prefixes = ["CAID",
                                    "HGVS",
                                    "MYVARIANT_HG38"]
# other options include ExAC.id or gnomAD.id, but we're not using those yet
curie_to_post_param_lookup = {
    "CAID": "id",
    "HGVS": "hgvs",
    "MYVARIANT_HG38": "MyVariantInfo_hg38.id"
}


class ClinGenService(object):

    logger = LoggingUtil.init_logging(__name__,
                                      logging.INFO,
                                      log_file_path=LoggingUtil.get_logging_path())

    def __init__(self):
        self.url = 'https://reg.genome.network/'
        self.synon_fields_param = 'fields=none+@id+externalRecords.dbSNP+externalRecords.ClinVarVariations+externalRecords.MyVariantInfo_hg38+genomicAlleles-genomicAlleles.referenceSequence'

    #
    # Important note: Provide a list of variant curies with the same prefix (ie. all HGVS or all CAID but not mixed)
    # That prefix must be one from the list: batchable_variant_curie_prefixes
    #
    def get_batch_of_synonyms(self, variant_curie_list: list):
        """
        Given a list of variant curies, return a corresponding list of sets of equivalent identifiers.

        ClinGenService.batchable_curie_prefixes provides a list of valid curie prefixes.

        :param variant_curie_list: a list of variant curies (with the same prefix)
        :return: a list of sets of equivalent identifiers - one for each variant curie supplied
        """
        if not variant_curie_list:
            return []

        curie_prefix = Text.get_curie(variant_curie_list[0])
        if curie_prefix not in batchable_variant_curie_prefixes:
            raise NotImplementedError(f'ClinGenService not able to support batches of {curie_prefix}!')

        variant_format_param = curie_to_post_param_lookup[curie_prefix]

        variant_id_list = [Text.un_curie(variant_curie) for variant_curie in variant_curie_list]

        separator = '\n'
        results_list = []
        batches = ceil(len(variant_id_list) / 2000)
        for i in range(batches):
            variant_subset = variant_id_list[i * 2000:i * 2000 + 2000]
            variant_pseudo_file = separator.join(variant_subset)
            query_url = f'{self.url}alleles?file={variant_format_param}&{self.synon_fields_param}'
            all_alleles_json = self.query_service(query_url, data=variant_pseudo_file)
            for index, allele_json in enumerate(all_alleles_json):
                results_list.append(self.parse_allele_json_for_synonyms(allele_json))
        return results_list

    def get_synonyms_by_other_id(self, variant_curie: str):
        variant_id_no_curie = Text.un_curie(variant_curie)
        variant_curie_prefix = Text.get_curie(variant_curie)
        if variant_curie_prefix.startswith('DBSNP'):
            possible_allele_preference = variant_id_no_curie.split("-")
            if len(possible_allele_preference) > 1:
                allele_preference = possible_allele_preference[1]
                variant_id_no_curie = possible_allele_preference[0]
            else:
                allele_preference = None
            return self.get_synonyms_by_parameter_matching('dbSNP.rs', variant_id_no_curie, allele_preference=allele_preference)

        elif variant_curie_prefix.startswith('CLINVARVARIANT'):
            return self.get_synonyms_by_parameter_matching('ClinVar.variationId', variant_id_no_curie)

        else:
            self.logger.debug(f'ClinGenService not able to support curie prefix {variant_curie_prefix}!')
            return None

    def get_synonyms_by_parameter_matching(self, url_param: str, url_param_value: str, allele_preference: str = None):
        synonym_sets = []
        query_url = f'{self.url}alleles?{url_param}={url_param_value}&{self.synon_fields_param}'
        query_json = self.query_service(query_url)
        for allele_json in query_json:
            synonym_sets.append(self.parse_allele_json_for_synonyms(allele_json))
        if allele_preference:
            filtered_syn_sets = []
            for syn_set in synonym_sets:
                robo_ids = Text.get_curies_by_prefix('ROBO_VARIANT', syn_set)
                if robo_ids:
                    actual_allele = robo_ids.pop().split('|')[-1]
                    if actual_allele == allele_preference:
                        filtered_syn_sets.append(syn_set)
            if filtered_syn_sets:
                return filtered_syn_sets

        return synonym_sets

    def parse_allele_json_for_synonyms(self, allele_json: dict):
        synonyms = set()
        try:
            variant_caid = allele_json['@id'].rsplit('/', 1)[1]
        except (KeyError, IndexError):
            return synonyms

        synonyms.add(f'CAID:{variant_caid}')

        if 'genomicAlleles' in allele_json:
            try:
                for genomic_allele in allele_json['genomicAlleles']:
                    for hgvs_id in genomic_allele['hgvs']:
                        synonyms.add(f'HGVS:{hgvs_id}')
                    if 'referenceGenome' in genomic_allele and genomic_allele['referenceGenome'] == 'GRCh38':
                        if 'chromosome' in genomic_allele:
                            sequence = genomic_allele['coordinates'][0]['allele']
                            chromosome = genomic_allele['chromosome']
                            start_position = genomic_allele['coordinates'][0]['start']
                            end_position = genomic_allele['coordinates'][0]['end']
                            robokop_variant_id = f'HG38|{chromosome}|{start_position}|{end_position}|{sequence}'
                            synonyms.add(f'ROBO_VARIANT:{robokop_variant_id}')

            except KeyError as e:
                error_message = f'parsing sequence variant synonym - genomicAlleles KeyError for {variant_caid}: {e}'
                self.logger.info(error_message)

        if 'externalRecords' in allele_json:
            if 'MyVariantInfo_hg19' in allele_json['externalRecords']:
                for myvar_json in allele_json['externalRecords']['MyVariantInfo_hg19']:
                    myvariant_id = myvar_json['id']
                    synonyms.add(f'MYVARIANT_HG19:{myvariant_id}')

            if 'MyVariantInfo_hg38' in allele_json['externalRecords']:
                for myvar_json in allele_json['externalRecords']['MyVariantInfo_hg38']:
                    myvariant_id = myvar_json['id']
                    synonyms.add(f'MYVARIANT_HG38:{myvariant_id}')

            if 'dbSNP' in allele_json['externalRecords']:
                for dbsnp_json in allele_json['externalRecords']['dbSNP']:
                    variant_rsid = dbsnp_json['rs']
                    synonyms.add(f'DBSNP:rs{variant_rsid}')

            if 'ClinVarVariations' in allele_json['externalRecords']:
                for clinvar_json in allele_json['externalRecords']['ClinVarVariations']:
                    clinvar_id = clinvar_json['variationId']
                    synonyms.add(f'CLINVARVARIANT:{clinvar_id}')

        return synonyms

    def get_variants_by_region(self, reference_sequence_label, center_position, region_size):
        flanking_size = int(region_size / 2)
        begin = center_position - flanking_size
        if begin < 0: begin = 0
        end = center_position + flanking_size
        query_url_main = f'{self.url}alleles?refseq={reference_sequence_label}&begin={begin}&end={end}&fields=none+@id&limit=2000&skip='
        counter = 0
        return_results = []
        while True:
            query_url = f'{query_url_main}{counter}'
            query_results = self.query_service(query_url)
            if query_results:
                for allele_json in query_results:
                    if '@id' in allele_json:
                        id_split = allele_json['@id'].rsplit('/', 1)
                        if (len(id_split) > 1) and ('CA' in id_split[1]):
                            variant_caid = id_split[1]
                            variant_node = SimpleNode(id=f'CAID:{variant_caid}', type=node_types.SEQUENCE_VARIANT)
                            return_results.append(variant_node)
                counter += 2000
            else:
                break
        return return_results

    def query_service(self, query_url, data=None, retries=1):
        if retries > 3:
            error_message = f'ClinGen Service tried but did not work for url: {query_url}'
            self.logger.error(error_message)
            raise requests.RequestException(error_message)
        try:
            if data:
                query_response = requests.post(query_url, data=data)
            else:
                query_response = requests.get(query_url)
            response_status_code = query_response.status_code
            if response_status_code != 200:
                if response_status_code == 404:
                    return {}
                else:
                    # this is probably a messed up request / our fault, skip retries
                    error_json = query_response.json()
                    cg_error_type = error_json["errorType"]
                    cg_error_description = error_json["description"]
                    cg_error_message = error_json["message"] if "message" in error_json else ""
                    error_message = f'ClinGen returned a 400 (bad request) calling ({query_url}):'
                    error_message += f'{cg_error_type} - {cg_error_description} - {cg_error_message}'
                    self.logger.error(error_message)
                    return self.query_service(query_url, data, 4)

            else:
                query_json = query_response.json()
                return query_json

        except requests.exceptions.RequestException as re:
            self.logger.error(f'Clingen service caught a request exception ({re}) on attempt number {retries}..')
            return self.query_service(query_url, data, retries + 1)

