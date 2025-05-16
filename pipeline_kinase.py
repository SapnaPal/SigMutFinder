import pandas as pd
import numpy as np
import glob
import os
import os.path
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import io
import gzip
import shutil
from scipy.stats import binomtest
from statsmodels.sandbox.stats.multicomp import multipletests
import gseapy as gp
import argparse
 
parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument("-maf_file", required=True, help="increase verbosity")

args = parser.parse_args()
config = vars(args)
print(config)

# Input file
maf = pd.read_csv(f'{args.maf_file}', sep='\t', encoding='windows-1252', low_memory=False)
maf = maf[['UNIPROT_ISOFORM', 'Hugo_Symbol', 'Variant_Classification', 'Amino_acids',
             'Protein_position', 'SIFT', 'IMPACT', 'hotspot', 'Tumor_Sample_Barcode']]

maf = maf.loc[maf['Variant_Classification'] == 'Missense_Mutation']

g = maf['Hugo_Symbol'].tolist()

# Pan-cancer Gene frequency

df = pd.read_csv('DATA/Pan_cancer_frequnecy_of_Kinase_genes.csv')
inp_gene_freq = df[df['Hugo_Symbol'].str.contains('|'.join(g), case=False)]
inp_gene_freq.to_csv('output_pan_cancer_kinase_gene_frequency.csv')


# Cancer-type gene frequency

df = pd.read_csv('DATA/Cancer_type_frequnecy_of_Kinase_genes.csv', index_col=0)
gene_data = df.columns.tolist()[1:]

input_gene = set(g)
data_gene = set(gene_data)

input_cancer_type_kianse_frequency = df[list(input_gene.intersection(data_gene))]

input_cancer_type_kianse_frequency.to_csv('output_cancer_type_kinase_gene_frequency.csv')

# Cancer-type mutation frequency

inp_kinase_gene = list(input_gene.intersection(data_gene))
maf['Mutation'] =  ' ' + maf['Amino_acids'].astype(str).str.split('/').str[0] + maf['Protein_position'].astype(str).str.split('/').str[0] + maf['Amino_acids'].astype(str).str.split('/').str[1]
maf1 = maf[maf['Hugo_Symbol'].str.contains('|'.join(inp_kinase_gene), case=False)]

inp_kinase_mut = [i.strip() for i in maf1['Mutation'].tolist()]
inp_kinase_tumor = input_cancer_type_kianse_frequency.T.columns.tolist()

file = []
for t in inp_kinase_tumor:
    tu = pd.read_csv(f'DATA/kinase_frequency/{t}_frequnecy_of_Kinase_genes.csv', index_col=0)
    d = tu[tu['Gene'].str.contains('|'.join(inp_kinase_gene), case=False) & tu['Mut'].str.contains('|'.join(inp_kinase_mut), case=False)]
    if len(d) > 0:
        file.append(d)

mut_freq = pd.concat(file)
mut_freq.to_csv('output_cancer_type_kinase_mutation_frequency.csv')

# Pathway enrichment 

df = pd.read_csv('DATA/kegg_enrichment_results.tsv', sep='\t', index_col=0)

path_enr = df[df['Genes'].str.contains('|'.join(inp_kinase_gene), case=False)]
path_enr.to_csv('output_enriched_pathway_having_kinase_gene.csv')