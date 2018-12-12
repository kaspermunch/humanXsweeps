

import sys
import mygene
import pandas
import numpy

def get_gene_coords(gene_list):
	mg = mygene.MyGeneInfo()

	records = list()
	for gene in mg.querymany(gene_list, scopes='symbol,accession,entrezgene,ensemblgene', 
		species='human', fields='symbol,genomic_pos_hg19'):

		if 'genomic_pos_hg19' not in gene:
			print('skipping', gene['query'], gene)
			continue

		coord = gene['genomic_pos_hg19']
		if not isinstance(coord, list):
			coord = [coord]
		for c in coord:
			records.append((gene['query'], gene['symbol'], c['chr'], int(c['start']), int(c['end'])))

	return pandas.DataFrame.from_records(records, columns=['query', 'symbol', 'chrom', 'start', 'end'])


som_file_name = '/Users/kmt/GenomeDK/simons/faststorage/people/kmt/data/chalmel_genes.csv'

hdf_file_name = '/Users/kmt/GenomeDK/simons/faststorage/people/kmt/results/chalmel_genes.hdf'

chalmel_genes = pandas.read_csv(som_file_name, sep=';').loc[lambda df: df.GeneName != '-']
chalmel_genes.loc[chalmel_genes.Pattern == '-'] = numpy.nan
chalmel_genes.Pattern = chalmel_genes.Pattern.astype(float)
chalmel_genes = chalmel_genes.loc[:, 'ID':'TF']


gene_list = chalmel_genes.GeneName.tolist()
coords = get_gene_coords(gene_list)

df = chalmel_genes.merge(coords, left_on=['GeneName'], right_on=['query'], how='left')

df.to_hdf(hdf_file_name, 'df', format='table', mode='w')

