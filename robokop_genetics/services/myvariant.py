from robokop_genetics.services.hgnc import HGNCService
from robokop_genetics.simple_graph_components import SimpleNode, SimpleEdge
from robokop_genetics.util import Text, LoggingUtil
from robokop_genetics import node_types

import requests
import logging
import time


class MyVariantService(object):

    def __init__(self, hgnc_service=None):
        log_file_path = LoggingUtil.get_logging_path()
        self.logger = LoggingUtil.init_logging(__name__,
                                               logging.INFO,
                                               log_file_path=log_file_path)
        self.url = "http://myvariant.info/v1/"
        self.effects_ignore_list = ['intergenic_region', 'sequence_feature']
        # we'll switch to this when they do
        #self.url_fields = 'snpeff.ann.effect,snpeff.ann.feature_type,snpeff.ann.gene_id'
        self.url_fields = 'snpeff.ann.effect,snpeff.ann.feature_type,snpeff.ann.genename'

        self.hgnc_service = hgnc_service if hgnc_service else HGNCService(log_file_path)

    def batch_sequence_variant_to_gene(self, variant_dict):
        if len(variant_dict) <= 1000:
            annotation_dictionary = {}
            post_params = {'fields': self.url_fields, 'ids': '', 'assembly': 'hg38'}
            id_lookup = {}
            for variant_id, variant_synonyms in variant_dict.items():
                # default to empty result for invalid or missing IDs
                annotation_dictionary[variant_id] = []
                # we could support hg19 as well, but calls need to be all one or the other
                # for now we only do hg38
                myvariant_curies = Text.get_curies_by_prefix('MYVARIANT_HG38', variant_synonyms)
                if not myvariant_curies:
                    self.logger.info(f'No MYVARIANT_HG38 synonym found for: {variant_id}')
                else:
                    for myvar_curie in myvariant_curies:
                        myvar_id = Text.un_curie(myvar_curie)
                        post_params['ids'] += f'{myvar_id},'
                        id_lookup[myvar_id] = variant_id

            if not post_params['ids']:
                self.logger.warning('batch_sequence_variant_to_gene called but all nodes provided had no MyVariant IDs')
                return annotation_dictionary

            # remove that extra comma
            post_params['ids'] = post_params['ids'][:-1]
            query_url = f'{self.url}variant'
            query_response = requests.post(query_url, data=post_params)
            if query_response.status_code == 200:
                query_json = query_response.json()
                for annotation_json in query_json:
                    try:
                        myvar_id = annotation_json['_id']
                        myvar_curie = f'MYVARIANT_HG38:{myvar_id}'
                        variant_id = id_lookup[myvar_id]
                        results = self.process_annotation(variant_id, annotation_json, myvar_curie)
                        if results:
                            annotation_dictionary[variant_id].extend(results)
                    except KeyError as e:
                        self.logger.warning(f'MyVariant batch call failed for annotation: {annotation_json["query"]} ({e})')
            else:
                self.logger.error(f'MyVariant batch non-200 response: ({query_response.status_code}) ids: ({post_params["ids"]})')
            return annotation_dictionary
        else:
            raise Exception('Error: More than 1000 variants not supported for MyVariant batch call.')
            return None

    def sequence_variant_to_gene(self, variant_id: str, variant_synonyms: set):
        return_results = []
        myvariant_ids = Text.get_curies_by_prefix('MYVARIANT_HG38', variant_synonyms)
        myvariant_assembly = 'hg38'
        # if we needed hg19
        #if not myvariant_ids:
        #    myvariant_ids = Text.get_curies_by_prefix('MYVARIANT_HG19', variant_synonyms)
        #    myvariant_assembly = 'hg19'
        if not myvariant_ids:
            self.logger.warning(f'No MyVariant ID found for {variant_id}, sequence_variant_to_gene failed.')
        else:
            for curie_myvariant_id in myvariant_ids:
                myvariant_id = Text.un_curie(curie_myvariant_id)
                query_url = f'{self.url}variant/{myvariant_id}?assembly={myvariant_assembly}&fields=snpeff'
                query_response = requests.get(query_url)
                if query_response.status_code == 200:
                    query_json = query_response.json()
                    return_results.extend(self.process_annotation(variant_id, query_json, curie_myvariant_id))
                else:
                    self.logger.error(f'MyVariant returned a non-200 response: {query_response.status_code})')

        return return_results

    def process_annotation(self, variant_id, annotation_json, curie_id):
        results = []
        try:
            if 'snpeff' in annotation_json:
                annotations = annotation_json['snpeff']['ann']
                # sometimes this is a list and sometimes a single instance
                if not isinstance(annotations, list):
                    annotations = [annotations]

                for annotation in annotations:
                    # for now we only take transcript feature type annotations
                    if annotation['feature_type'] != 'transcript':
                        continue

                    # TODO: this assumes the gene_id is also a HGNC symbol and not an ID
                    # for now we look and find real ID if we can
                    # in the future we'd like to do the following
                    # gene_identifier = f'HGNC:{annotation["gene_id"]}'
                    gene_symbol = annotation['genename']
                    gene_id = self.hgnc_service.get_gene_id_from_symbol(gene_symbol)
                    if gene_id is None:
                        # if we can't find a real id, just skip it
                        self.logger.info(f'Could not find real ID for gene symbol: {gene_symbol}')
                        continue

                    gene_node = SimpleNode(id=gene_id, name=gene_symbol, type=node_types.GENE)

                    props = {}
                    #do we want this?
                    #if 'putative_impact' in annotation:
                    #    props['putative_impact'] = annotation['putative_impact']

                    effects_list = annotation['effect'].split('&')
                    for effect in effects_list:
                        if effect in self.effects_ignore_list:
                            continue

                        predicate_id = f'SNPEFF:{effect}'
                        predicate_label = effect

                        edge = SimpleEdge(source_id=variant_id,
                                          target_id=gene_node.id,
                                          provided_by='myvariant.sequence_variant_to_gene',
                                          input_id=curie_id,
                                          predicate_id=predicate_id,
                                          predicate_label=predicate_label,
                                          ctime=time.time(),
                                          properties=props)
                        results.append((edge, gene_node))
            
            else:
                self.logger.error(f'No snpeff annotation found for variant {variant_id}')

        except KeyError as e:
            self.logger.error(f'Myvariant annotation error:{e}')
                
        return results
