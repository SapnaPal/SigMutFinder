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


# Pan-cancer mutation information
df1 = pd.read_csv('DATA/pan_cancer_total_mutation.csv', index_col=0)


#Input file
maf = pd.read_csv(f'{args.maf_file}', sep='\t', encoding='windows-1252', low_memory=False)
maf = maf[['UNIPROT_ISOFORM', 'Hugo_Symbol', 'Variant_Classification', 'Amino_acids',
             'Protein_position', 'SIFT']]

maf = maf.loc[maf['Variant_Classification'] == 'Missense_Mutation']

g = maf['Hugo_Symbol'].tolist()

# Significant PPIs having input file genes
sg = pd.read_csv('DATA/pan_cancer_ppi_significance.csv', index_col=0)
sg1 = sg.loc[sg['fdr_value'] <= 0.05]

sg1['p1_int_mutations_uniq'] = sg1['p1_int_mutations_name'].apply(lambda x: list(set(x[1:-1].split(','))))
sg1['p2_int_mutations_uniq'] = sg1['p2_int_mutations_name'].apply(lambda x: list(set(x[1:-1].split(','))))

sg1['p1_int_mutations_uniq_total'] = sg1['p1_int_mutations_name'].apply(lambda x: len(list(set(x[1:-1].split(',')))))
sg1['p2_int_mutations_uniq_total'] = sg1['p2_int_mutations_name'].apply(lambda x: len(list(set(x[1:-1].split(',')))))
sg1 = sg1[sg1['P1_gene'].str.contains('|'.join(g), case=False) | sg1['P2_gene'].str.contains('|'.join(g), case=False)]

# Input genes involved in significant PPI interactions
p1 = set(sg1['P1_gene'].to_list())
p2 = set(sg1['P2_gene'].to_list())
p3  = set(list(p1.union(p2)))
p4 = set(g)

# Input file with significant PPI genes
maf1 = maf[maf['Hugo_Symbol'].str.contains('|'.join(list(p3.intersection(p4))), case=False)]

# Pan-cancer mutation information having input genes
a = df1[df1['Hugo_Symbol'].str.contains('|'.join(list(p3.intersection(p4))), case=False)]

# Significant PPIs involving input genes with mutations, p-value, and fdr-value information

final_inp_map = []
for q in maf1['Hugo_Symbol'].tolist():
    inp =  maf1[maf1['Hugo_Symbol'] == q]
    #inp = maf1.iloc[q]
    gene = inp['Hugo_Symbol'].tolist()[0]
    sg2 = sg1[sg1['P1_gene'].str.contains(gene, case=False) | sg1['P2_gene'].str.contains(gene, case=False)]
        
    gene1 = []
    protein_position = []
    amino_acids = []



    if len(sg2) >0:
        for i in range(0, len(sg2)):
            data = sg2.iloc[i]
        
            inp_mut = inp['Protein_position'].tolist()[0]
        
            #p1
            data_mut1 = data['p1_int_mutations_uniq']
            for x in data_mut1:
                if inp_mut in x:
                    print(i)
                    print(f'P1 {inp_mut}')
                    gene1.append(inp['Hugo_Symbol'].tolist()[0])
                    protein_position.append(inp_mut)
                    amino_acids.append(inp['Amino_acids'].tolist()[0])
                

            #p2
            data_mut2 = data['p2_int_mutations_uniq']
            for x in data_mut2:
                if inp_mut in x:
                    print(i)
                    print(f'P2 {inp_mut}')
                    gene1.append(inp['Hugo_Symbol'].tolist()[0])
                    protein_position.append(inp_mut)
                    amino_acids.append(inp['Amino_acids'].tolist()[0])

    if len(gene1) >0:
        
        a1 = a.loc[a['Hugo_Symbol'] == gene1[0]]
        a2 = a1.loc[a1['Protein_position'] == protein_position[0]]
        a2 = a1.loc[a1['Amino_acids'] == amino_acids[0]]
        a2['cancer'].tolist()

        ppi_name = []
        pvalue = []
        fdrvalue = []

        for i in range(0 , len(sg2)):
            t = sg2.iloc[i]
            ppi = f"{t['P1_gene']}_{t['P2_gene']}"
            p_value = f"{t['p_value']}"
            fdr_value = f"{t['fdr_value']}"
            ppi_name.append(ppi)
            pvalue.append(p_value)
            fdrvalue.append(fdr_value)

        inp['cancer'] = [list(set(a2['cancer'].tolist()))]
        inp['ppi_name'] = [ppi_name]
        inp['pvalue'] = [pvalue]
        inp['fdrvalue'] = [fdrvalue]
        final_inp_map.append(inp)
        
result = pd.concat(final_inp_map)
result.to_csv('output.csv')