from ftplib import FTP, error_proto, all_errors as all_ftp_errors
from io import BytesIO
from json import loads
import logging
import time
from robokop_genetics.util import LoggingUtil


def pull_via_ftp(ftpsite, ftpdir, ftpfile):
    ftp = FTP(ftpsite)
    ftp.login()
    ftp.cwd(ftpdir)
    with BytesIO() as data:
        ftp.retrbinary(f'RETR {ftpfile}', data.write)
        binary = data.getvalue()
    ftp.quit()
    return binary


class HGNCService(object):

    def __init__(self, log_file_path=None):
        if log_file_path is None:
            self.logging_on = False
        else:
            self.logging_on = True
            self.logger = LoggingUtil.init_logging(__name__,
                                                   logging.INFO,
                                                   logFilePath=log_file_path)
        self.hgnc_symbol_to_curie = None

    def get_gene_id_from_symbol(self, gene_symbol: str):
        if self.hgnc_symbol_to_curie is None:
            self.init_symbol_lookup()
        if gene_symbol in self.hgnc_symbol_to_curie:
            return self.hgnc_symbol_to_curie[gene_symbol]
        else:
            if self.logging_on:
                self.logger.info(f'HGNCService could not find ID for gene symbol: {gene_symbol}')
            return None

    def init_symbol_lookup(self):
        if self.logging_on:
            self.logger.info(f'Preparing HGNC Symbol look up.')
        self.hgnc_symbol_to_curie = {}

        data = None
        num_tries = 0
        while num_tries < 5 and data is None:
            try:
                if self.logging_on:
                    self.logger.info(f'Pulling HGNC set by ftp.')
                data = pull_via_ftp('ftp.ebi.ac.uk', '/pub/databases/genenames/new/json', 'hgnc_complete_set.json')
            except all_ftp_errors:
                num_tries += 1
                time.sleep(2)
                if self.logging_on:
                    self.logger.info(f'FTP attempt failed. Trying again ({num_tries} times).')
        if data is None:
            if self.logging_on:
                self.logger.info(f'HGNC Symbol look up failed.!')
            return

        hgnc_json = loads(data.decode())
        for hgnc_item in hgnc_json['response']['docs']:
            hgnc_symbol = hgnc_item['symbol']
            if hgnc_symbol not in self.hgnc_symbol_to_curie:
                self.hgnc_symbol_to_curie[hgnc_symbol] = hgnc_item['hgnc_id']

        if self.logging_on:
            self.logger.info(f'HGNC Symbol look up ready.!')