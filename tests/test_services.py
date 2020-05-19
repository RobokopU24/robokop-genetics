import pytest

from robokop_genetics.genetics_services import *


@pytest.fixture()
def genetics_services():
    # can use this for troubleshooting
    #log_file_path = './'
    log_file_path = None
    return GeneticsServices(log_file_path)


def test_gene_symbol_to_id(genetics_services):
    gene_id = genetics_services.get_gene_id_from_symbol('ASS1')
    assert gene_id == 'HGNC:758'

    gene_id = genetics_services.get_gene_id_from_symbol('DMD')
    assert gene_id == 'HGNC:2928'

    gene_id = genetics_services.get_gene_id_from_symbol('BRCA1')
    assert gene_id == 'HGNC:1100'

    gene_id = genetics_services.get_gene_id_from_symbol('THISISAFAKEGENE')
    assert gene_id is None


def test_myvariant(genetics_services):

    relations = genetics_services.get_variant_to_gene(MYVARIANT, 'MYVARIANT_HG19:chr7:g.55241707G>T', {'MYVARIANT_HG19:chr7:g.55241707G>T'})
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:3236' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'missense_variant' in plabels
    assert 'downstream_gene_variant' in plabels
    pids = [edge.predicate_id for edge, node in relations]
    assert 'SNPEFF:missense_variant' in pids
    assert 'SNPEFF:downstream_gene_variant' in pids

    relations = genetics_services.get_variant_to_gene(MYVARIANT, 'MYVARIANT_HG38:chr11:g.68032291C>G', {'MYVARIANT_HG38:chr11:g.68032291C>G'})
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:7715' in identifiers
    assert 'HGNC:41796' in identifiers
    assert 'HGNC:410' in identifiers

    relations = genetics_services.get_variant_to_gene(MYVARIANT, 'MYVARIANT_HG19:chr17:g.56283533T>A', {'MYVARIANT_HG19:chr17:g.56283533T>A'})
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:7121' in identifiers
    assert 'HGNC:3423' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'splice_acceptor_variant' in plabels
    assert 'downstream_gene_variant' in plabels
    pids = [edge.predicate_id for edge, node in relations]
    assert 'SNPEFF:splice_acceptor_variant' in pids
    assert 'SNPEFF:downstream_gene_variant' in pids


def test_batch_myvariant(genetics_services):

    variant_dict = dict()
    variant_dict['MYVARIANT_HG38:chr11:g.68032291C>G'] = {'MYVARIANT_HG38:chr11:g.68032291C>G'}
    variant_dict['MYVARIANT_HG38:chrX:g.32389644G>A'] = {'MYVARIANT_HG38:chrX:g.32389644G>A'}
    variant_dict['MYVARIANT_HG38:chr17:g.7674894G>A'] = {'MYVARIANT_HG38:chr17:g.7674894G>A'}
    variant_dict['MYVARIANT_HG38:chr9:g.130489423A>G'] = {'MYVARIANT_HG38:chr9:g.130489423A>G'}

    batch_results = genetics_services.batch_variant_to_gene(MYVARIANT, variant_dict)

    relations = batch_results['MYVARIANT_HG38:chr11:g.68032291C>G']
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:7715' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'missense_variant' in plabels

    relations = batch_results['MYVARIANT_HG38:chrX:g.32389644G>A']
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:2928' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'stop_gained' in plabels

    relations = batch_results['MYVARIANT_HG38:chr17:g.7674894G>A']
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:11998' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'stop_gained' in plabels

    relations = batch_results['MYVARIANT_HG38:chr9:g.130489423A>G']
    identifiers = [node.id for edge, node in relations]
    assert 'HGNC:758' in identifiers
    plabels = [edge.predicate_label for edge, node in relations]
    assert 'missense_variant' in plabels


def test_ensembl(genetics_services):
    # using hg38
    robokop_variant_id = f'ROBO_VARIANT:HG38|17|58206171|58206172|A'
    relations = genetics_services.get_variant_to_gene(ENSEMBL, 'CAID:CA279509', {robokop_variant_id})
    identifiers = [node.id for edge, node in relations]
    assert 'ENSEMBL:ENSG00000011143' in identifiers
    assert 'ENSEMBL:ENSG00000121053' in identifiers
    assert 'ENSEMBL:ENSG00000167419' in identifiers
    assert len(identifiers) > 20

    robokop_variant_id = f'ROBO_VARIANT:HG38|1|69092|69093|C'
    relations = genetics_services.get_variant_to_gene(ENSEMBL, 'CAID:CA16728208', {robokop_variant_id})
    identifiers = [node.id for edge, node in relations]
    assert 'ENSEMBL:ENSG00000186092' in identifiers
    assert 'ENSEMBL:ENSG00000240361' in identifiers

    robokop_variant_id = f'ROBO_VARIANT:HG38|1|11796321|11796322|A'
    relations = genetics_services.get_variant_to_gene(ENSEMBL, 'CAID:CA170990', {robokop_variant_id})
    identifiers = [node.id for edge, node in relations]
    assert 'ENSEMBL:ENSG00000177000' in identifiers
    assert 'ENSEMBL:ENSG00000011021' in identifiers

    robokop_variant_id = f'ROBO_VARIANT:HG38|X|32389643|32389644|A'
    relations = genetics_services.get_variant_to_gene(ENSEMBL, 'CAID:CA267021', {robokop_variant_id})
    identifiers = [node.id for edge, node in relations]
    assert 'ENSEMBL:ENSG00000198947' in identifiers
