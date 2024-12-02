from robokop_genetics import node_types
from robokop_genetics.simple_graph_components import SimpleNode, SimpleEdge
from robokop_genetics.util import Text, LoggingUtil
from collections import namedtuple
import logging
import sqlite3
import os
import requests

###
# !!!! README - this uses ensembl gene information to annotate variants with nearby genes,
# but it hasn't been used in a long time. It's been replaced by running SNPEFF in ORION.
###

EnsemblGene = namedtuple('EnsemblGene', ['ensembl_id', 'ensembl_name', 'chromosome', 'start_position', 'end_position', 'gene_biotype', 'description'])


class EnsemblService(object):

    logger = LoggingUtil.init_logging(__name__,
                                      logging.INFO,
                                      log_file_path=LoggingUtil.get_logging_path())
    
    def __init__(self, temp_dir: str=None):

        self.upstream_gene_predicate_id = 'SNPEFF:upstream_gene_variant'
        self.upstream_gene_predicate_label = 'upstream_gene_variant'

        self.downstream_gene_predicate_id = 'SNPEFF:downstream_gene_variant'
        self.downstream_gene_predicate_label = 'downstream_gene_variant'

        # This can cause issues if python doesn't have write access to the directory.
        # If so, an appropriate directory needs to be specified. See the README.
        temp_dir = temp_dir if temp_dir else '.'
        self.gene_db_path = os.path.join(temp_dir, 'genes.sqlite3')

        self.gene_db_successfully_created = False
        self.persistent_conn = None
        self.all_gene_annotations = None

        # we assume the order of attributes from this url -
        # if we change this we need to change the indexing in create_genes_db below
        self.ensembl_genes_url = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
                                    <!DOCTYPE Query>
                                    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                                        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                                            <Attribute name = "ensembl_gene_id" />
                                            <Attribute name = "gene_biotype" />
                                            <Attribute name = "external_gene_name" />
                                            <Attribute name = "start_position" />
                                            <Attribute name = "end_position" />
                                            <Attribute name = "description" />
                                            <Attribute name = "chromosome_name" />
                                        </Dataset>
                                    </Query>"""

        self.check_if_already_done_sql = "SELECT name FROM sqlite_master WHERE type='table' AND name='genes';"
        
        self.genes_table_sql = """CREATE TABLE IF NOT EXISTS genes (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        ensembl_id text, 
        gene_name text, 
        chromosome text,
        start_pos INTEGER, 
        end_pos INTEGER, 
        gene_type text,
        description text);"""

        self.genes_table_ensembl_id_index_sql = "CREATE UNIQUE INDEX ensembl_ids on genes(ensembl_id);"
        self.genes_table_composite_index_sql = "CREATE INDEX gene_composite on genes(chromosome, start_pos, end_pos, ensembl_id);"

        self.gene_entry_sql = """INSERT INTO genes 
        (ensembl_id, gene_name, chromosome, start_pos, end_pos, gene_type, description) 
        VALUES (?,?,?,?,?,?,?);"""
        
        self.gene_range_select_sql = """SELECT ensembl_id, gene_name, start_pos, end_pos
        FROM genes WHERE chromosome = ? AND ((? >= start_pos AND ? <= end_pos)
        OR (? >= start_pos AND ? <= end_pos) OR (? <= start_pos AND ? >= end_pos));"""

        self.gene_ensembl_id_select_sql = "SELECT * FROM genes WHERE ensembl_id = ?"

    def create_or_connect_to_genes_db(self):
        
        if not self.gene_db_successfully_created:
            self.create_genes_db()
        
        if not self.persistent_conn:
            self.persistent_conn = sqlite3.connect(self.gene_db_path, detect_types=sqlite3.PARSE_DECLTYPES)
            self.persistent_conn.row_factory = sqlite3.Row

        return self.persistent_conn

    def create_genes_db(self):
        try:
            db_conn = sqlite3.connect(self.gene_db_path)
            db_cursor = db_conn.cursor()

            db_cursor.execute(self.check_if_already_done_sql)
            if db_cursor.fetchone() is not None:
                # TODO we assume this means it's complete - it could have been created but not filled out
                self.gene_db_successfully_created = True
                return True

            # ensembl_genes are of type EnsemblGene namedtuple - that has to match the SQL parameters in gene_entry_sql
            ensembl_genes = self.retrieve_all_genes()
           
            with db_conn:
                db_conn.execute(self.genes_table_sql)
                db_conn.executemany(self.gene_entry_sql, ensembl_genes)
                db_conn.execute(self.genes_table_ensembl_id_index_sql)
                db_conn.execute(self.genes_table_composite_index_sql)

            db_conn.close()
            self.logger.info(f'Ensembl created a gene database with {len(ensembl_genes)} entries!')
            self.gene_db_successfully_created = True
            return True

        except sqlite3.Error as e:
            self.logger.error(f'Ensembl had a database error: {e}')

    def retrieve_all_genes(self):
        genes_response = requests.get(self.ensembl_genes_url)
        if genes_response.status_code == 200:
            genes_data = genes_response.text.splitlines()
            if len(genes_data) > 1:
                ensembl_genes = []
                for gene_line in genes_data[1:]:
                    # gene_info is a EnsemblGene 
                    gene_info = self.parse_biomart_gene_data(gene_line)
                    if gene_info:
                        ensembl_genes.append(gene_info)
                return ensembl_genes
            else:
                self.logger.error(f'Ensembl biomart genes call didnt find any matches! Thats not right!')
                return False
        else:
            self.logger.error(f'Ensembl non-200 response from biomart genes call: {genes_response.status_code})')
            return None

    def parse_biomart_gene_data(self, gene_line):
        gene_data = gene_line.split('\t')
        try:
            # we assume the index order is predictable from ensembl_genes_url above - if we change that we need to change this
            ensembl_id = gene_data[0]
            gene_type = gene_data[1]
            ensembl_name = gene_data[2]
            start_position = int(gene_data[3])
            end_position = int(gene_data[4])
            description = gene_data[5]
            chromosome = gene_data[6]
            return EnsemblGene(ensembl_id, ensembl_name, chromosome, start_position, end_position, gene_type, description)

        except (IndexError, ValueError) as e:
            self.logger.error(f'Ensembl biomart genes call had an issue with one line: {e})')

        return None

    def sequence_variant_to_gene(self, variant_id: str, variant_synonyms: set):
        
        flanking_region_size = 500000
        results = []

        found_valid_robokop_key = False
        robokop_ids = Text.get_curies_by_prefix('ROBO_VARIANT', variant_synonyms)
        if not robokop_ids:
            self.logger.debug(f'ensembl: robokop variant key not found for variant: {variant_id}')
            return results
        else:
            try:
                for robokop_key in robokop_ids:
                    robokop_data = Text.un_curie(robokop_key).split('|')
                    reference_genome = robokop_data[0]
                    if reference_genome == 'HG38':
                        found_valid_robokop_key = True
                        robokop_key_used = robokop_key
                    else:
                        continue
                    chromosome = robokop_data[1]
                    start_position = int(robokop_data[2])
                    end_position = int(robokop_data[3])
            except IndexError as e:
                self.logger.error(f'ensembl: robokop variant key not set properly for variant: {variant_id} - {robokop_ids[0]}')
                return results

        if not found_valid_robokop_key:
            self.logger.debug(f'ensembl: latest robokop variant key not found for variant: {variant_id}')
            return results

        flanking_min = start_position - flanking_region_size
        if flanking_min < 0:
            flanking_min = 0
        flanking_max = end_position + flanking_region_size

        db_conn = self.create_or_connect_to_genes_db()
        db_cursor = db_conn.cursor()

        #logger.info(f'looking for genes overlapping {flanking_min}-{flanking_max}')

        db_cursor.execute(self.gene_range_select_sql, (chromosome, flanking_min, flanking_min, flanking_max, flanking_max, flanking_min, flanking_max))

        genes_in_region = db_cursor.fetchall()
        for gene_id_text, gene_name, gene_start, gene_end in genes_in_region:
            #cast this to make neo4j happy
            gene_id = str(gene_id_text)
            #logger.info(f'Found matching gene: {gene_id},{gene_start},{gene_end}')
            gene_node = SimpleNode(id=f'ENSEMBL:{gene_id}', name=f'{gene_name}', type=node_types.GENE)
            if start_position < gene_start:
                distance = gene_start - start_position
                predicate_id = self.upstream_gene_predicate_id
                predicate_label = self.upstream_gene_predicate_label
            elif end_position > gene_end:
                distance = end_position - gene_end
                predicate_id = self.downstream_gene_predicate_id
                predicate_label = self.downstream_gene_predicate_label
            else:
                distance = 0
                predicate_id = self.downstream_gene_predicate_id
                predicate_label = self.downstream_gene_predicate_label

            props = {'distance': distance}
            edge = SimpleEdge(source_id=variant_id,
                              target_id=gene_node.id,
                              provided_by='ensembl.sequence_variant_to_gene',
                              input_id=robokop_key_used,
                              predicate_id=predicate_id,
                              predicate_label=predicate_label,
                              ctime=1,
                              properties=props)
            results.append((edge, gene_node))

        self.logger.debug(f'ensembl sequence_variant_to_gene found {len(results)} results for {variant_id}')

        return results

    def get_all_ensembl_gene_annotations(self):
        if self.all_gene_annotations:
            return self.all_gene_annotations
        else:
            ensembl_genes = self.retrieve_all_genes()
            all_gene_annotations = {}
            for gene in ensembl_genes:
                all_gene_annotations[gene.ensembl_id] = {
                    'name': gene.ensembl_name,
                    'chromosome': gene.chromosome,
                    'start_position': gene.start_position,
                    'end_position': gene.end_position,
                    'gene_biotype': gene.gene_biotype,
                    'description': gene.description
                }
            self.all_gene_annotations = all_gene_annotations
            return all_gene_annotations

    def get_ensembl_gene_annotations(self, ensembl_id):
        if not self.all_gene_annotations:
            self.get_all_ensembl_gene_annotations()

        if ensembl_id in self.all_gene_annotations:
            return self.all_gene_annotations[ensembl_id]
        else:
            return {'ensembl_error': f'{ensembl_id} not found.'}
