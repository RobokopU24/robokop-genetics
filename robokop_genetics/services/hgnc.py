import json
import logging
import time
import requests
from robokop_genetics.util import LoggingUtil


class HGNCService(object):

    logger = LoggingUtil.init_logging(__name__,
                                      logging.INFO,
                                      log_file_path=LoggingUtil.get_logging_path())

    def __init__(self):
        self.hgnc_symbol_to_curie = None

    def get_gene_id_from_symbol(self, gene_symbol: str):
        if self.hgnc_symbol_to_curie is None:
            self.init_symbol_lookup()
        if gene_symbol in self.hgnc_symbol_to_curie:
            return self.hgnc_symbol_to_curie[gene_symbol]
        else:
            self.logger.info(f'HGNCService could not find ID for gene symbol: {gene_symbol}')
            return None

    def init_symbol_lookup(self):
        self.logger.debug(f'Preparing HGNC Symbol look up.')
        self.hgnc_symbol_to_curie = {}

        hgnc_json = None
        num_tries = 0
        while num_tries < 5 and hgnc_json is None:
            try:
                self.logger.debug(f'Pulling HGNC data.')
                hgnc_response = requests.get("https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json")
                try:
                    hgnc_json = hgnc_response.json()
                except json.JSONDecodeError:
                    self.logger.error(f'HGNC download json parsing error.')
            except requests.exceptions.RequestException:
                num_tries += 1
                time.sleep(2)
                self.logger.warning(f'HGNC download attempt failed. Trying again ({num_tries} times).')
        if hgnc_json is None:
            self.logger.error(f'HGNC Symbol look up failed.!')
            return

        for hgnc_item in hgnc_json['response']['docs']:
            hgnc_symbol = hgnc_item['symbol']
            if hgnc_symbol not in self.hgnc_symbol_to_curie:
                self.hgnc_symbol_to_curie[hgnc_symbol] = hgnc_item['hgnc_id']
        self.logger.debug(f'HGNC Symbol look up ready.!')
