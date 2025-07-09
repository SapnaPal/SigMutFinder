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


parser.add_argument("-maf_file_path", required=True, help="increase verbosity")

args = parser.parse_args()
config = vars(args)
print(config)




#total sample information
def fcount(path):
    count1 = 0
    for root, dirs, files in os.walk(path):
            count1 += len(dirs)

    return count1

path = f'{args.maf_file_path}'
tumor_type = []
tumor_samples = []
for x in glob.glob(f'{path}*'):
    tumor = x.split("/")[-1]
    tumor_type.append(tumor)
    _, _, files = next(os.walk(f'{path}{tumor}/'))
    tumor_samples.append(len(files))
    
print(tumor_samples)  #total samples download pan-cancer


#remove duplicated samples by barcode matching

cancer_maf_files = []

for t in tumor_type:    
    all_filesl = glob.glob(os.path.join(f'{args.maf_file_path}{t}/', "*.maf"))
    cancer = []
    files = []
    barcode = []
    for a in all_filesl:
        f = pd.read_csv(a, sep='\t', encoding='windows-1252', low_memory=False)
        if len(f) >0:
            cancer.append(t)
            files.append(a)
            f1 = f["Tumor_Sample_Barcode"].loc[0]
            barcode.append(f1)
    data = {'cancer': cancer, 'file': files, 'barcode_full': barcode}
    tumor_data = pd.DataFrame(data)
    b = tumor_data 
    b['barcode'] = b['barcode_full'].str[:12]
    b = b.drop_duplicates(subset='barcode', keep="last")
    cancer_maf_files.append(b)
pan_cancer_data = pd.concat(cancer_maf_files)


#total samples pan-cancer

sample_counts = []
maf_data_list = []
for t in tumor_type:
    all_files = pan_cancer_data.loc[pan_cancer_data['cancer'] == t]['file'].tolist()
    df = pd.concat((pd.read_csv(f, sep='\t', encoding='windows-1252', low_memory=False) for f in all_files)
                   , ignore_index=True)
    df['Cancer_type'] = t
    maf_data_list.append(df)
    sample_counts.append(df["Tumor_Sample_Barcode"].nunique())
    
print(sum(sample_counts))


#Tumor sample table for plotting

data = {'Cancer type': tumor_type, 'Total sample': sample_counts}
tumor_data = pd.DataFrame(data)
tumor_data.head(5)

tumor_data.to_csv('cancer_type_sample_info.csv')
print('cancer_type_sample_info saved')

# Prepare pan-cancer consolidated mutation file for kinase specific analysis
maf_data = pd.concat(maf_data_list, names=["Cancer_Type", "Original_Index"])
maf_data_final = maf_data[maf_data['Variant_Classification'].str.contains('|'.join(['Frame_Shift_Del', 'Frame_Shift_Ins', 
                                                                   'In_Frame_Del', 'In_Frame_Ins', 
                                                                   'Missense_Mutation', 'Nonsense_Mutation']))]

#total samples pan-cancer

mut_counts = []
cancer_type = []
for t in tumor_type:
    mut_counts.append(len(maf_data_final[maf_data_final["Cancer_type"] == t]))
    cancer_type.append(t)
    
sum(sample_counts)


#Tumor sample table for plotting

data = {'Cancer type': cancer_type, 'Total mutations': mut_counts}
tumor_data = pd.DataFrame(data)
tumor_data.to_csv('cancer_type_mutation_info.csv')
print('cancer_type_mutation_info saved')

# Total number of unique samples across all cancer types
total_samples = maf_data_final["Tumor_Sample_Barcode"].nunique()
print('Total number of unique samples across all cancer types')
print(total_samples)


#SIFT scores
df1 = maf_data_final
df1 = df1['SIFT'].dropna()
df1['sift_score'] = maf_data_final['SIFT'].str.split('(').str[1].str[:-1]
df1['sift_score'] = df1['sift_score'].astype(float)
print(f'Mutatins with SIFT scores {len(df1)}')

x = df1.index.tolist()
y = df1['sift_score'].tolist()

sns.set_theme(rc={'figure.figsize':(35,20)}, style='white')
sns.set(context='notebook', style='white', palette='deep', font='sans-serif', font_scale=4, color_codes=False, rc=None)

sns.kdeplot(data=y, color="turquoise", label="EGFR-inactive-Wild", linewidth=2, fill=True)
plt.xlabel("SIFT-score")
plt.ylabel("Density of SIFT scores for mutations")
plt.savefig('cohort_sift_score_of_pan_cancer_mutations.png', bbox_inches='tight')


#kinase gene info
kinase_info = pd.read_csv('../uniprotkb_keyword_KW_0418_AND_reviewed_2024_04_03.tsv', sep='\t')
#kinase_info

kinase = []
for i in range(0, 636):
    x = kinase_info['Gene Names'][i]
    x = x.replace(' ', ',')
    ki = x.split(',')[0]
    kinase.append(ki)   
gene_list = kinase 

# Pan-cancer consolidated mutation file specific to kinase genes
maf_data_kinase = maf_data_final[maf_data_final["Hugo_Symbol"].isin(gene_list)]
maf_data_kinase_save = maf_data_kinase[['UNIPROT_ISOFORM', 'Hugo_Symbol', 'Variant_Classification', 'Amino_acids',
             'Protein_position', 'SIFT', 'Cancer_type']]
maf_data_kinase_save.to_csv('cohort_Mutations_in_kinases.csv')
print('mutations_present_in_kinases file saved')

#Total number of tumor samples with specific Kinase mutation
mutation_frequency = maf_data_kinase.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()
mutation_frequency.to_csv('cohort_Total_number_of_tumor _samples_with_specific_Kinase_mutation.csv')

#Pan_cancer_frequnecy_of_Kinase_genes
mutation_frequency_percentage = (mutation_frequency / total_samples) * 100
mutation_frequency_percentage = mutation_frequency_percentage.sort_values(ascending=False)
mutation_frequency_percentage.to_csv("cohort_Pan_cancer_frequnecy_of_Kinase_genes.csv")

#Cancer_type_frequnecy_of_Kinase_genes
cancer_type_frequencies = maf_data_kinase.groupby(["Cancer_type", "Hugo_Symbol"])["Tumor_Sample_Barcode"].nunique().unstack(fill_value=0)
cancer_type_frequencies_percentage = (cancer_type_frequencies.T / cancer_type_frequencies.sum(axis=1)) * 100
cancer_type_frequencies_percentage = cancer_type_frequencies_percentage.T
cancer_type_frequencies_percentage.to_csv("cohort_Cancer_type_frequnecy_of_Kinase_genes.csv")

frequency = []
tumor  = []
Gene = []

# Top frequent kinase genes (Cancer-type specific)
for t in tumor_type:
    fre = round(cancer_type_frequencies_percentage.loc[[t], :].T.sort_values(by=t, ascending=False).iloc[0].tolist()[0], 2)
    frequency.append(fre)
    Gene.append(cancer_type_frequencies_percentage.loc[[t], :].T.sort_values(by=t, ascending=False).index.tolist()[0])
    tumor.append(t)
data = {'tumor': tumor, 'Gene': Gene, 'frequency': frequency}
freq_data = pd.DataFrame(data)
freq_data['kinase1'] = freq_data['tumor'].astype(str) + '\n' + '(' + freq_data['frequency'].astype(str) + '%' + ')'
freq_data = freq_data.sort_values(['Gene'])
freq_data.to_csv('cohort_top_frequent_kinase_genes_Cancer_type_specific.csv')


combined_maf = maf_data_kinase
combined_maf['Combined'] =  ' ' + combined_maf['Amino_acids'].astype(str).str.split('/').str[0] + combined_maf['Protein_position'].astype(str).str.split('/').str[0] + combined_maf['Amino_acids'].astype(str).str.split('/').str[1]

freq_mut = []
for t in tumor_type:
    cancer_df = combined_maf[combined_maf['Cancer_type'].str.contains(t, na=False)]
    cancer_df['Mutation_Key'] = cancer_df[['Hugo_Symbol', 'Variant_Classification', 
                                       'Combined']].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    # Count the occurrences of each mutation
    mutation_counts = cancer_df['Mutation_Key'].value_counts().reset_index()
    mutation_counts.columns = ['Mutation', 'Count']

    # Sort by frequency (descending)
    mutation_counts = mutation_counts.sort_values(by='Count', ascending=False)
    total_mutations = mutation_counts["Count"].sum()
    mutation_counts["Frequency (%)"] = (mutation_counts["Count"] / total_mutations) * 100
    mutation_counts["tumor"] = t
    mutation_counts["Mut"] = mutation_counts['Mutation'].astype(str).str.split(' ').str[1]
    mutation_counts["Gene"] = mutation_counts['Mutation'].astype(str).str.split('_').str[0]
    mutation_counts.to_csv(f"cohort_{t}_frequnecy_of_mutations_in_Kinase_genes.csv")   #save kinase mutation freq for each cancer type
    freq_mut.append(mutation_counts.head(1))
    
can_freq_mut = pd.concat(freq_mut)
can_freq_mut['kinase1'] = can_freq_mut['tumor'].astype(str) + '\n' + can_freq_mut['Mut'].astype(str) + '\n' +  '(' + can_freq_mut['Frequency (%)'].round(2).astype(str) + '%' + ')'
can_freq_mut['freq'] = can_freq_mut['Frequency (%)'].round(2)
df = can_freq_mut.sort_values(['Gene'])
df.to_csv(f"cohort_top_frequent_mutation_of_Kinase_gene_cancer_types.csv")

