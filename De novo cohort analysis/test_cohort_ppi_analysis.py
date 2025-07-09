import pandas as pd
import numpy as np
import glob
import os
import os.path
#import matplotlib.pyplot as plt
#import seaborn as sns
import gzip
import io
import gzip
import shutil
from scipy.stats import binomtest
from statsmodels.sandbox.stats.multicomp import multipletests
import warnings
import argparse


# Suppress FutureWarning messages
warnings.simplefilter(action='ignore', category=FutureWarning)


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


eclare = pd.read_csv('../H_sapiens_interfacesHQ.txt', sep='\t')
filtered_df = eclare[eclare['P1_IRES'] != '[]'] 
filtered_ppi_data = filtered_df[filtered_df['P2_IRES'] != '[]']

p1 = set(filtered_ppi_data['P1'].to_list())
p2 = set(filtered_ppi_data['P2'].to_list())
p3  = list(p1.union(p2))
dfp3 = pd.DataFrame({'gene': p3})
#dfp3.to_csv('all_ppi_proteins.csv', index=False)
#len(set(p3))


protein = pd.read_csv('../mutation_mapping/idmapping_reviewed_true_AND_model_organ_2024_04_13.tsv', sep='\t') #uniprot gene info
pkinase = []
for i in range(0, 8689):
    x = protein['Gene Names'][i]

    x = str(x).replace(' ', ',')
    ki = x.split(',')[0]
    pkinase.append(ki)
    #print(ki)
protein['gene'] = pkinase


ppi_id = protein['From'].to_list()

p1 = filtered_ppi_data['P1'].to_list()
p1 = list(set(p1))
p2 = filtered_ppi_data['P2'].to_list()
p2 = list(set(p2))
removep1 = list(set(p1) - set(ppi_id))
removep2 = list(set(p2) - set(ppi_id))
#print(len(set(p1) - set(ppi_id)))
#print(len(set(p2) - set(ppi_id)))


filtered_ppi_data = filtered_ppi_data[~filtered_ppi_data['P1'].str.contains('|'.join(removep1), case=False)]
filtered_ppi_data = filtered_ppi_data[~filtered_ppi_data['P2'].str.contains('|'.join(removep2), case=False)]
filtered_ppi_data = filtered_ppi_data.reset_index()


gene1 = []
length1 = []
for i in range(0,61532):
    f = filtered_ppi_data.iloc[i:i+1]
    x = f['P1'][i]
    y = protein[protein['From'].str.contains(x, na=False)]
    geney = y['gene'].to_list()
    gene1.append(geney[0])
    lengthy = y['Length'].to_list()
    length1.append(lengthy[0])


filtered_ppi_data['P1_gene'] = gene1
filtered_ppi_data['P1_length'] = length1


gene2 = []
length2 = []
for i in range(0,61532):
    f = filtered_ppi_data.iloc[i:i+1]
    x = f['P2'][i]
    y = protein[protein['From'].str.contains(x, na=False)]

    geney = y['gene'].to_list()
    #print(x)
    gene2.append(geney[0])

    lengthy = y['Length'].to_list()
    length2.append(lengthy[0])

filtered_ppi_data['P2_gene'] = gene2
filtered_ppi_data['P2_length'] = length2



p1_intereface_info = []
p2_intereface_info = []

p1_non_intereface_info = []
p2_non_intereface_info = []


p1_interface_length = []
p2_interface_length = []

p1_non_interface_length = []
p2_non_interface_length = []


for h in range(0,61532):
    prot1 = []
    x = filtered_ppi_data['P1_IRES'][h]
    df1 = filtered_ppi_data['P1_length'][h]
    df2 = df1.astype(int)

    x = x[1:-1]

    for i in x.split(','):
        if  "-" in i:
            e = i.split("-")
            f = int(e[0])
            g = int(e[1])
            h = range(f, g+1)
            for i1 in h:
                prot1.append(int(i1))
        else:
            prot1.append(int(i))
    p1_intereface_info.append(prot1)
    p1_interface_length.append(len(prot1))


    l1 = []
    for j in range(1, df2+1):
        l1.append(j)
    non_interface = []
    for k in l1:
        if k not in prot1:
            non_interface.append(k)
    p1_non_intereface_info.append(non_interface)
    p1_non_interface_length.append(len(non_interface))

for h in range(0,61532):
    prot2 = []
    x = filtered_ppi_data['P2_IRES'][h]
    df1 = filtered_ppi_data['P2_length'][h]
    df2 = df1.astype(int)

    x = x[1:-1]

    for i in x.split(','):
        if  "-" in i:
            e = i.split("-")
            f = int(e[0])
            g = int(e[1])
            h = range(f, g+1)
            for i1 in h:
                prot2.append(int(i1))
        else:
            prot2.append(int(i))
    p2_intereface_info.append(prot2)
    p2_interface_length.append(len(prot2))

    l1 = []
    for j in range(1, df2+1):
        l1.append(j)
    non_interface = []
    for k in l1:
        if k not in prot2:
            non_interface.append(k)
    p2_non_intereface_info.append(non_interface)
    p2_non_interface_length.append(len(non_interface))
    
    
filtered_ppi_data['P1_Intfc'] = p1_intereface_info
filtered_ppi_data['P2_intfc'] = p2_intereface_info
filtered_ppi_data['P1_non_intfc'] = p1_non_intereface_info
filtered_ppi_data['P2_non_intfc'] = p2_non_intereface_info
filtered_ppi_data['P1_intfc_length'] = p1_interface_length
filtered_ppi_data['P2_intfc_length'] = p2_interface_length
filtered_ppi_data['P1_non_intfc_length'] = p1_non_interface_length
filtered_ppi_data['P2_non_intfc_length'] = p2_non_interface_length
ppi_kinase = filtered_ppi_data.drop(['P1_IRES', 'P2_IRES'], axis=1)

all_files = pan_cancer_data['file'].tolist()
df = pd.concat((pd.read_csv(f, sep='\t', encoding='windows-1252', low_memory=False) for f in all_files)
                   , ignore_index=True)
cancer_type = pan_cancer_data['cancer'].tolist()

all_files1 = []
for f,c in zip(all_files, cancer_type):
    f1 = pd.read_csv(f, sep='\t', encoding='windows-1252', low_memory=False)
    f1['cancer'] = c
    all_files1.append(f1)
    
df = pd.concat(all_files1)

name = [] #names of the kinases
total = [] #total mutations found for each kinase in cancer type LUAD
uniq_mut = []
int_mut = []
non_int_mut = []
df1 = df[['cancer', 'UNIPROT_ISOFORM', 'Hugo_Symbol', 'Variant_Classification', 'Amino_acids', 
             'Protein_position', 'SIFT', 'IMPACT', 'hotspot', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode']]

m = ppi_kinase

p1_inM = []
p1_non_inM = []
p2_inM = []
p2_non_intM = []
p1_m_total = []
p2_m_total = []


p1_int_mutation_name = []
p1_int_mutation_can_name = []
p2_int_mutation_name = []
p2_int_mutation_can_name = []

p1_nonint_mutation_name = []
p1_nonint_mutation_can_name = []
p2_nonint_mutation_name = []
p2_nonint_mutation_can_name = []

for s in range(0,61532):
    
    m = ppi_kinase
    m = m.iloc[s:s+1]
    P1_gene = m['P1_gene'].to_list()
    P1_l = m['P1_length'].to_list()
    P1_in_res = m['P1_Intfc'].to_list()
    P1_nonint_res = m['P1_non_intfc'].to_list()

    P2_gene = m['P2_gene'].to_list()
    P2_l = m['P2_length'].to_list()
    P2_in_res = m['P2_intfc'].to_list()
    P2_nonint_res = m['P2_non_intfc'].to_list()

    P1_in_m = []   #correct 0
    P1_nonint_m = []  #correct 2
    P2_in_m = []    #correct 0
    P2_nonint_m = []   #correct 1

    P1_total_mut = []   #correct 2
    P2_total_mut = []   #correct 1
    P1_in_m_name = []
    P1_in_m_can_name = []
    P2_in_m_name = []
    P2_in_m_can_name = []
    
    P1_nonint_m_name = []
    P1_nonint_m_can_name = []
    P2_nonint_m_name = []
    P2_nonint_m_can_name = []
    

    for i1,j1, i2, j2 in zip(P1_gene, P1_l, P2_gene, P2_l):


        p1_df2 = df1.loc[df1['Hugo_Symbol'] == i1]
        p1_df3 = p1_df2[p1_df2['SIFT'].str.contains('deleterious', na=False)]
        p1_ppb = p1_df3[p1_df3['Protein_position'].str.contains(f'/{j1}', na=False)]
        
       
        p2_df2 = df1.loc[df1['Hugo_Symbol'] == i2]
        p2_df3 = p2_df2[p2_df2['SIFT'].str.contains('deleterious', na=False)]
        p2_ppb = p2_df3[p2_df3['Protein_position'].str.contains(f'/{j2}', na=False)]

        if len(p1_ppb) and len(p2_ppb) >0:
            
            #print(len(p1_ppb), len(p2_ppb))
            p1_pp = p1_ppb['Protein_position'].to_list()
            p1_aa = p1_ppb['Amino_acids'].to_list()
           
            P1_total_mut.append(len(p1_ppb))
            p1_int_m = []
            p1_non_int_m = []
            p1_in_m_name = []
            p1_in_m_can = []
            p1_nonint_m_name = []
            p1_nonint_m_can = []
            
            
            for i in p1_pp:
                l = [int(i.split('/')[0])]
                
                
               
                
                for k in P1_in_res[0]:
                    if k in l:
                        
                        p1_int_m.append(k)
                        
                        p1_p = p1_ppb[p1_ppb['Protein_position'].str.contains(f'{k}/', na=False)]
                        aa = p1_p['Amino_acids'].to_list()
                        p_p = p1_p['Protein_position'].to_list()
                        
                        P1_name =  aa[0] + ' ' + p_p[0]

                        p1_in_m_name.append(P1_name)

                        p1_cancer = p1_p['cancer'].to_list()
                        for i_can in p1_cancer:
                            #l_can = i_can[5:7]
                            p1_in_m_can.append(i_can)
              
                        
                        
                for k in P1_nonint_res[0]:
                    if k in l:
                        
                        p1_non_int_m.append(k)
                        #aa = p1_ppb['Amino_acids'].to_list()
                       # p_p = p1_ppb['Protein_position'].to_list()
                       
                        p1_p = p1_ppb[p1_ppb['Protein_position'].str.contains(f'{k}/', na=False)]
                        aa = p1_p['Amino_acids'].to_list()
                        p_p = p1_p['Protein_position'].to_list()
                        
                        P1_name =  aa[0] + ' ' + p_p[0]
                        
                        p1_nonint_m_name.append(P1_name)
                        
                        p1_cancer = p1_p['cancer'].to_list()
                        for i_can in p1_cancer:
                            #l_can = i_can[5:7]
                            p1_nonint_m_can.append(i_can)
                        
                        
                        
            P1_in_m.append(len(p1_int_m))
            P1_nonint_m.append(len(p1_non_int_m))
           
            
            p2_pp = p2_ppb['Protein_position'].to_list()
            p2_aa = p2_ppb['Amino_acids'].to_list()
           
            P2_total_mut.append(len(p2_ppb))
            p2_int_m = []
            p2_non_int_m = []
            p2_in_m_name = []
            p2_in_m_can = []
            p2_nonint_m_name = []
            p2_nonint_m_can = []
            
            for i in p2_pp:
                l = [int(i.split('/')[0])]
                
                
                
                for k in P2_in_res[0]:
                    if k in l:
                        p2_int_m.append(k)
                        
                        p2_p = p2_ppb[p2_ppb['Protein_position'].str.contains(f'{k}/', na=False)]
                        aa2 = p2_p['Amino_acids'].to_list()
                        p_p2 = p2_p['Protein_position'].to_list()
                        
                        P2_name =  aa2[0] + ' ' + p_p2[0]

                        p2_in_m_name.append(P2_name)

                        p2_cancer = p2_p['cancer'].to_list()
                        for i_can2 in p2_cancer:
                            #l_can2 = i_can2[5:7]
                            p2_in_m_can.append(i_can2)



                for k in P2_nonint_res[0]:
                    if k in l:
                        p2_non_int_m.append(k)
                        
                        p2_p = p2_ppb[p2_ppb['Protein_position'].str.contains(f'{k}/', na=False)]
                        aa2 = p2_p['Amino_acids'].to_list()
                        p_p2 = p2_p['Protein_position'].to_list()
                        
                        P2_name =  aa2[0] + ' ' + p_p2[0]
                        #print(P2_name)
                        
                        p2_nonint_m_name.append(P2_name)
                        
                        p2_cancer = p2_p['cancer'].to_list()
                        for i_can2 in p2_cancer:
                            #l_can2 = i_can2[5:7]
                            p2_nonint_m_can.append(i_can2)
                        
                        
                        
                        
            P2_in_m.append(len(p2_int_m))
            P2_nonint_m.append(len(p2_non_int_m))
            
            
            
            P1_in_m_name.append(p1_in_m_name)
            P1_in_m_can_name.append(p1_in_m_can)
            P2_in_m_name.append(p2_in_m_name)
            P2_in_m_can_name.append(p2_in_m_can)
            
            P1_nonint_m_name.append(p1_nonint_m_name)
            P1_nonint_m_can_name.append(p1_nonint_m_can)
            P2_nonint_m_name.append(p2_nonint_m_name)
            P2_nonint_m_can_name.append(p2_nonint_m_can)
            
            #print(p2_nonint_m_name)
            #print(p2_nonint_m_can)
            
        else:
            #print('Nan')
            P1_in_m.append('N')
            P1_nonint_m.append('N')
            P2_in_m.append('N')
            P2_nonint_m.append('N')
            
            P1_total_mut.append('N')
            P2_total_mut.append('N')

            P1_in_m_name.append('N')
            P1_in_m_can_name.append('N')
            P2_in_m_name.append('N')
            P2_in_m_can_name.append('N')
            
            P1_nonint_m_name.append('N')
            P1_nonint_m_can_name.append('N')
            P2_nonint_m_name.append('N')
            P2_nonint_m_can_name.append('N')
    
    print(P1_nonint_m_name[0])
    p1_inM.append(P1_in_m[0])
    p1_non_inM.append(P1_nonint_m[0])
    p2_inM.append(P2_in_m[0])
    p2_non_intM.append(P2_nonint_m[0])
    p1_m_total.append(P1_total_mut[0])
    p2_m_total.append(P2_total_mut[0])

    p1_int_mutation_name.append(list(P1_in_m_name[0]))
    p1_int_mutation_can_name.append(list(P1_in_m_can_name[0]))
    p2_int_mutation_name.append(list(P2_in_m_name[0]))
    p2_int_mutation_can_name.append(list(P2_in_m_can_name[0]))
    
    p1_nonint_mutation_name.append(list(P1_nonint_m_name[0]))
    p1_nonint_mutation_can_name.append(list(P1_nonint_m_can_name[0]))
    p2_nonint_mutation_name.append(list(P2_nonint_m_name[0]))
    p2_nonint_mutation_can_name.append(list(P2_nonint_m_can_name[0]))
    
test = ppi_kinase
test['P1_int_m'] = p1_inM
test['P1_non_int_m'] = p1_non_inM
test['P2_int_m'] = p2_inM
test['P2_non_int_m'] = p2_non_intM


test['p1_int_mutations_name'] = p1_int_mutation_name
test['p1_int_mutation_tumor_name'] = p1_int_mutation_can_name
test['p2_int_mutations_name'] = p2_int_mutation_name
test['p2_int_mutation_tumor_name'] = p2_int_mutation_can_name


test['p1_nonint_mutations_name'] = p1_nonint_mutation_name
test['p1_nonint_mutation_tumor_name'] = p1_nonint_mutation_can_name
test['p2_nonint_mutations_name'] = p2_nonint_mutation_name
test['p2_nonint_mutation_tumor_name'] = p2_nonint_mutation_can_name




test['P1_total_mutations'] = p1_m_total
test['P2_total_mutations'] = p2_m_total


source = test[~test.P1_int_m.str.contains("N", na=False)]
print('Final file prepared for.')
binomial1 = []
binomial2 = []
p_value = []

source = source.loc[source['P1_total_mutations'] != 0]
source = source.loc[source['P2_total_mutations'] != 0]
print(len(source))


for s in range(0,len(source)):
    interface1 = source.iloc[s]
    x = interface1['P1_int_m']
    n = interface1['P1_total_mutations']
    p = interface1['P1_intfc_length']/interface1['P1_length']
    p = binomtest(x, n=n, p=p, alternative='greater')
    q = p.pvalue
    binomial1.append(q)

    interface2 = source.iloc[s]
    a = interface2['P2_int_m']
    b = interface2['P2_total_mutations']
    c = interface1['P2_intfc_length']/interface2['P2_length']
    d = binomtest(a, n=b, p=c, alternative='greater')
    e = d.pvalue
    binomial2.append(e)

    t1 = q*e
    p_value.append(t1)

    
# Create a list of the adjusted p-values
p_adjusted = multipletests(p_value, alpha=0.05, method='bonferroni')

source['p_value'] = p_value
source['fdr_value'] = p_adjusted[1]
f = source.sort_values(by=['fdr_value'], ascending=True)
f1 = f.reset_index()


f1.to_csv('cohort_pan_cancer_ppi_significance.csv')
print('Pan cancer Interface significance saved.')