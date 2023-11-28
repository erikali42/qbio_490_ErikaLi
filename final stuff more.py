#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
os.chdir('/Users/erikali/Desktop/QBIO/qbio_490_ErikaLi/analysis_data')


# In[2]:


import cptac
cptac.download(dataset="Endometrial")
endometrial = cptac.Endometrial()


# In[3]:


protein_data = endometrial.get_proteomics()
protein_data.columns = protein_data.columns.get_level_values(0) 
protein_data


# In[38]:


clinical_data = endometrial.get_clinical()
clinical_data


# In[5]:


clinical_data['Num_full_term_pregnancies'].value_counts()


# In[6]:


rna_data = endometrial.get_transcriptomics()
rna_data


# In[46]:


merged_data = pd.merge(protein_data, clinical_data, on='Patient_ID')
merged_data


# In[47]:


columns_to_remove = ['Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',
    'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN',
    'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM', 'tumor_Stage-Pathological',
    'FIGO_stage', 'LVSI', 'BMI','Sample_ID', 'Sample_Tumor_Normal', 'Proteomics_Tumor_Normal', 'Country',
    'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive'
]

merged_data.drop(columns=columns_to_remove, inplace=True)


# In[48]:


merged_data


# In[49]:


merged_data['Num_full_term_pregnancies'].value_counts()


# In[79]:


import numpy as np
import pandas as pd

# Handle 'Num_full_term_pregnancies' column
merged_data['Num_full_term_pregnancies'].replace('4 or more', '4', inplace=True)
merged_data['Num_full_term_pregnancies'] = merged_data['Num_full_term_pregnancies'].replace({'Unknown': np.nan}).astype(float)

# Create masks for low and high pregnancy groups
lowpregnancy_mask = merged_data['Num_full_term_pregnancies'].isin([1.0, 2.0])
highpregnancy_mask = merged_data['Num_full_term_pregnancies'].isin([3.0, 4.0])

# Filter data based on masks
lowpregnancy_proteins = merged_data[lowpregnancy_mask]
highpregnancy_proteins = merged_data[highpregnancy_mask]

# Calculate mean expression levels
mean_lowpregnancy = lowpregnancy_proteins.mean()
mean_highpregnancy = highpregnancy_proteins.mean()

# Calculate absolute difference in expression
expression_difference = np.abs(mean_lowpregnancy - mean_highpregnancy)

# Sort proteins based on absolute expression difference
sorted_proteins = expression_difference.sort_values(ascending=False).dropna()

# Extract top 5 proteins with the largest expression differences
top_50_proteins = sorted_proteins.head(50)

print(top_50_proteins)


# In[80]:


from bioservices import KEGG

# Instantiate KEGG service
k = KEGG()

# List of genes
gene_list = [
    'SPINK6', 'KRTAP4-1', 'KRTAP4-6', 'XAGE1A', 'TCHHL1', 'PAGE2B', 'KRT28', 'KRTAP5-2', 'BPIFA1', 'PAGE5',
    'FAM25C', 'IL36RN', 'SPRR2G', 'KRT71', 'CADM3', 'KRTAP9-2', 'KRT25', 'S100A7A', 'KRT27', 'BSX', 'RBP2',
    'CPLX2', 'C4orf46', 'KRTAP3-3', 'LIN28A', 'SST', 'ABO', 'TPH1', 'DSG4','LGALS4', 'KRTDAP', 'CHST11', 'STMN4', 'TAF1L', 'PMCH', 'RBMXL3', 'IL36G', 'CRNN', 'H2AFB2', 'DSG3',
    'SERPINA11', 'LOC441155', 'GAGE5', 'IFI6', 'KRTAP4-9', 'KRT34', 'ITLN2', 'PTGES3L-AARSD1', 'PPAN-P2RY11'
]

# Retrieve KEGG pathway information for your genes
kegg_pathways = {}
for gene in gene_list:
    pathway = k.get_pathway_by_gene(gene, organism='hsa')
    kegg_pathways[gene] = pathway

# Display pathway information
for gene, pathway_info in kegg_pathways.items():
    print(f"Pathways for {gene}:")
    print(pathway_info)
    print("------------------------")


# In[81]:


from gprofiler import GProfiler

# Set your user agent string (replace 'Your_User_Agent' with an appropriate identifier)
user_agent = 'Your_User_Agent'

# Initialize g:Profiler with the user agent
gp = GProfiler(user_agent=user_agent)

# List of genes you want to perform enrichment analysis on
genes_of_interest = [
    'SPINK6', 'KRTAP4-1', 'KRTAP4-6', 'XAGE1A', 'TCHHL1', 'PAGE2B', 'KRT28', 'KRTAP5-2', 'BPIFA1', 'PAGE5',
    'FAM25C', 'IL36RN', 'SPRR2G', 'KRT71', 'CADM3', 'KRTAP9-2', 'KRT25', 'S100A7A', 'KRT27', 'BSX', 'RBP2',
    'CPLX2', 'C4orf46', 'KRTAP3-3', 'LIN28A', 'SST', 'ABO', 'TPH1', 'DSG4', 'Num_full_term_pregnancies',
    'LGALS4', 'KRTDAP', 'CHST11', 'STMN4', 'TAF1L', 'PMCH', 'RBMXL3', 'IL36G', 'CRNN', 'H2AFB2', 'DSG3',
    'SERPINA11', 'LOC441155', 'GAGE5', 'IFI6', 'KRTAP4-9', 'KRT34', 'ITLN2', 'PTGES3L-AARSD1', 'PPAN-P2RY11'
]

# Perform pathway enrichment analysis
enrichment_results = gp.gprofile(organism='hsapiens', query=genes_of_interest)

# Display pathway enrichment results
print(enrichment_results)


# In[83]:


import matplotlib.pyplot as plt

# Extracting relevant data for plotting
data = [
    [1, True, 3.5e-08, 216, 10, 7, 0.7, 0.032, 'REAC:R-HSA-6805567', 'rea', 1, 'Keratinization', 1, 'KRT71,SPRR2G,KRT28,SPINK6,KRTAP4-6,KRTAP4-1,KRTAP5-2'],
    [1, True, 7.21e-08, 224, 13, 7, 0.538, 0.031, 'GO:0031424', 'BP', 1, 'keratinization', 1, 'KRT71,SPRR2G,KRT28,SPINK6,KRTAP4-6,KRTAP4-1,KRTAP5-2'],
    [1, True, 5.46e-07, 299, 13, 7, 0.538, 0.023, 'GO:0030216', 'BP', 1, 'keratinocyte differentiation', 1, 'KRT71,SPRR2G,KRT28,SPINK6,KRTAP4-6,KRTAP4-1,KRTAP5-2'],
    [1, True, 1.74e-06, 353, 13, 7, 0.538, 0.02, 'GO:0009913', 'BP', 1, 'epidermal cell differentiation', 1, 'KRT71,SPRR2G,KRT28,SPINK6,KRTAP4-6,KRTAP4-1,KRTAP5-2'],
    [1, True, 5.09e-06, 412, 13, 7, 0.538, 0.017, 'GO:0043588', 'BP', 1, 'skin development', 1, 'KRT71,SPRR2G,KRT28,SPINK6,KRTAP4-6,KRTAP4-1,KRTAP5-2']
    

# Extracting p-values and processes
p_values = [d[2] for d in data]  # Negating p-values for better visualization
processes = [d[11] for d in data]

# Creating the plot
plt.figure(figsize=(10, 6))
plt.barh(processes, p_values, color='skyblue')
plt.xlabel('P-values')
plt.title('Functional Enrichment Analysis')
plt.gca().invert_yaxis()  # Invert y-axis to have highest values at the top

# Set y-axis ticks and labels to biological processes
plt.yticks(processes, processes)

# Set x-axis limits based on the range of p-values
plt.xlim(min(p_values), max(p_values))
plt.xticks(p_values, ['%.1e' % val for val in p_values])

plt.tight_layout()

# Show plot
plt.show()


# In[90]:


import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Your gProfiler output data (example subset)
data = [
    1, True, 2.29e-16, 216, 28, 15, 0.536, 0.069, 'REAC:R-HSA-6805567', 'rea', 1, 'Keratinization', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.97e-15, 224, 45, 15, 0.333, 0.067, 'GO:0031424', 'BP', 1, 'keratinization', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.06e-13, 299, 45, 15, 0.333, 0.05, 'GO:0030216', 'BP', 1, 'keratinocyte differentiation', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.61e-12, 353, 45, 15, 0.333, 0.042, 'GO:0009913', 'BP', 1, 'epidermal cell differentiation', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 7.42e-12, 459, 45, 16, 0.356, 0.035, 'GO:0008544', 'BP', 1, 'epidermis development', 1, 'RBP2,KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.51e-11, 412, 45, 15, 0.333, 0.036, 'GO:0043588', 'BP', 1, 'skin development', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 6.02e-10, 29, 15, 6, 0.4, 0.207, 'HPA:020040', 'hpa', 1, 'hair; cells in internal root sheath', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 6.02e-10, 29, 15, 6, 0.4, 0.207, 'HPA:020020', 'hpa', 1, 'hair; cells in cuticle', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020000', 'hpa', 1, 'hair', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020030', 'hpa', 1, 'hair; cells in external root sheath', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020010', 'hpa', 1, 'hair; cells in cortex/medulla', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3']


# Extracting relevant data for plotting
enrichment_ratios = [d[6] for d in data]
biological_processes = [d[11] for d in data]
fdr_values = [d[7] for d in data]

# Creating the plot
plt.figure(figsize=(10, 8))
bars = plt.barh(biological_processes, enrichment_ratios, color=cm.viridis(fdr_values))

# Adding labels and titles
plt.xlabel('Enrichment Ratio')
plt.ylabel('Biological Pathway')
plt.title('Functional Enrichment Analysis')

# Adding a colorbar separately
sm = cm.ScalarMappable(cmap=cm.viridis)
sm.set_array(fdr_values)
cbar = plt.colorbar(sm, label='FDR')

# Show plot
plt.tight_layout()
plt.show()


# In[92]:


import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Your gProfiler output data (example subset)
data = [
    1, True, 2.29e-16, 216, 28, 15, 0.536, 0.069, 'REAC:R-HSA-6805567', 'rea', 1, 'Keratinization', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.97e-15, 224, 45, 15, 0.333, 0.067, 'GO:0031424', 'BP', 1, 'keratinization', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.06e-13, 299, 45, 15, 0.333, 0.05, 'GO:0030216', 'BP', 1, 'keratinocyte differentiation', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.61e-12, 353, 45, 15, 0.333, 0.042, 'GO:0009913', 'BP', 1, 'epidermal cell differentiation', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 7.42e-12, 459, 45, 16, 0.356, 0.035, 'GO:0008544', 'BP', 1, 'epidermis development', 1, 'RBP2,KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 3.51e-11, 412, 45, 15, 0.333, 0.036, 'GO:0043588', 'BP', 1, 'skin development', 1, 'KRT34,DSG3,KRT71,SPRR2G,KRT27,KRT28,DSG4,SPINK6,KRTAP4-6,KRTAP4-1,KRT25,KRTAP5-2,KRTAP4-9,KRTAP3-3,KRTAP9-2'], [1, True, 6.02e-10, 29, 15, 6, 0.4, 0.207, 'HPA:020040', 'hpa', 1, 'hair; cells in internal root sheath', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 6.02e-10, 29, 15, 6, 0.4, 0.207, 'HPA:020020', 'hpa', 1, 'hair; cells in cuticle', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020000', 'hpa', 1, 'hair', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020030', 'hpa', 1, 'hair; cells in external root sheath', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'], [1, True, 7.52e-10, 30, 15, 6, 0.4, 0.2, 'HPA:020010', 'hpa', 1, 'hair; cells in cortex/medulla', 1, 'KRT34,KRT71,KRT27,KRT28,DSG4,KRTAP3-3'],[1, True, 0.0107, 136, 16, 4, 0.25, 0.029, 'KEGG:04915', 'keg', 1, 'Estrogen signaling pathway', 1, 'KRT34,KRT27,KRT28,KRT25']



# Extracting relevant data for plotting
enrichment_ratios = [d[6] for d in data]
biological_processes = [d[11] for d in data]
fdr_values = [d[7] for d in data]

# Creating the plot
plt.figure(figsize=(10, 8))
bars = plt.barh(biological_processes, enrichment_ratios, color=cm.viridis(fdr_values))

# Adding labels and titles
plt.xlabel('Enrichment Ratio')
plt.ylabel('Biological Pathway')
plt.title('Functional Enrichment Analysis')

# Adding a colorbar separately
sm = cm.ScalarMappable(cmap=cm.viridis)
sm.set_array(fdr_values)
cbar = plt.colorbar(sm, label='FDR')

# Show plot
plt.tight_layout()
plt.show()


# In[58]:


import seaborn as sns
import matplotlib.pyplot as plt

# Assuming you've filtered lowpregnancy_proteins and highpregnancy_proteins

# Normalize the data (if needed)
normalized_lowpregnancy = (lowpregnancy_proteins - lowpregnancy_proteins.mean()) / lowpregnancy_proteins.std()
normalized_highpregnancy = (highpregnancy_proteins - highpregnancy_proteins.mean()) / highpregnancy_proteins.std()

# Create separate heatmaps for low and high pregnancy groups
fig, axs = plt.subplots(1, 2, figsize=(18, 8))  # 1 row, 2 columns for separate heatmaps

# Heatmap for low pregnancy
sns.heatmap(normalized_lowpregnancy, cmap='viridis', ax=axs[0], annot=False)
axs[0].set_title('Low Pregnancy Protein Expression')
axs[0].set_xlabel('Proteins')
axs[0].set_ylabel('Samples')

# Heatmap for high pregnancy
sns.heatmap(normalized_highpregnancy, cmap='viridis', ax=axs[1], annot=False)
axs[1].set_title('High Pregnancy Protein Expression')
axs[1].set_xlabel('Proteins')
axs[1].set_ylabel('Samples')

plt.tight_layout()
plt.show()


# In[69]:


from scipy.stats import ttest_ind

# Calculate mean expression for low and high pregnancy
mean_low = lowpregnancy_proteins.mean()
mean_high = highpregnancy_proteins.mean()

# Avoid division by zero or NaNs
mean_low_nonzero = mean_low.replace(0, np.nan)
mean_high_nonzero = mean_high.replace(0, np.nan)

# Calculate fold change considering non-zero values for low pregnancy
fold_change_low = np.log2(mean_high_nonzero / mean_low_nonzero)

# Perform t-test for each protein for low pregnancy
t_test_results_low = []
for protein in lowpregnancy_proteins.columns:
    t_stat, p_val = ttest_ind(lowpregnancy_proteins[protein], highpregnancy_proteins[protein])
    t_test_results_low.append((protein, p_val))

# Extract p-values for low pregnancy
p_values_low = [result[1] for result in t_test_results_low]

# Repeat the same procedure for high pregnancy
# Calculate fold change considering non-zero values for high pregnancy
fold_change_high = np.log2(mean_low_nonzero / mean_high_nonzero)

# Perform t-test for each protein for high pregnancy
t_test_results_high = []
for protein in highpregnancy_proteins.columns:
    t_stat, p_val = ttest_ind(highpregnancy_proteins[protein], lowpregnancy_proteins[protein])
    t_test_results_high.append((protein, p_val))

# Extract p-values for high pregnancy
p_values_high = [result[1] for result in t_test_results_high]


# In[65]:


import matplotlib.pyplot as plt
import numpy as np

# Define alpha threshold for significance
alpha = 0.05

# Calculate the negative log of p-values
neg_log_p_values_low = -np.log10(p_values_low)
neg_log_p_values_high = -np.log10(p_values_high)

# Plotting the volcano plot for low pregnancy
plt.figure(figsize=(8, 6))
plt.scatter(fold_change_low, neg_log_p_values_low, color='blue', alpha=0.5)

# Highlight significantly different proteins (based on alpha threshold)
significant_low = np.where(np.array(p_values_low) < alpha)
plt.scatter(fold_change_low[significant_low], neg_log_p_values_low[significant_low], color='red', alpha=0.5)
plt.xlabel('Fold Change (log2)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot for Low Pregnancy')
plt.axhline(-np.log10(alpha), color='green', linestyle='--', label=f'p-value threshold: {alpha}')
plt.legend()
plt.show()

# Plotting the volcano plot for high pregnancy
plt.figure(figsize=(8, 6))
plt.scatter(fold_change_high, neg_log_p_values_high, color='blue', alpha=0.5)
# Highlight significantly different proteins (based on alpha threshold)
significant_high = np.where(np.array(p_values_high) < alpha)
plt.scatter(fold_change_high[significant_high], neg_log_p_values_high[significant_high], color='red', alpha=0.5)
plt.xlabel('Fold Change (log2)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot for High Pregnancy')
plt.axhline(-np.log10(alpha), color='green', linestyle='--', label=f'p-value threshold: {alpha}')
plt.legend()
plt.show()


# In[ ]:




