#!/usr/bin/env python
# coding: utf-8

# # Analysis of the Mapping Data
# 
# Input
# 
# -> Neurommsig_pd_full_genes.csv
# 
# - A list of Biological mechanisms (NeuroMMSig) along with genes in those mechanisms associated to PD is contained in 
# 
# -> RiskSNPs.rda
# 
# - A mapping from the above-mentioned genes to SNPs in the respective mechanisms is available in the RData object
# 
# # Goals / Statistical relevant questions
# Map Mechanisms to Gene to SNPs
# 
#  - How many Mechanism can actually be mapped to genes that reveal prs scores from the PheWAS catalog. How many 
#      Mechanisms are useful then?
#  - For how many genes can we calculate prs scores regarding the mapped snps?
#  - Side note: Whats the difference between comp and compDF?

# In[102]:


# imports


# In[1]:


import pyreadr
import math
import pandas as pd


# In[2]:


def prs(genes,csv,log=False,p_threshold=5.000000e-02):
    
    """ Function that takes the genes list + csv as input and filters the csv by the given genes and "Parkinson's disease" phenotype.
        Optional is a p-value threshold (default=0.05) that dismisses all entries larger than the set threshold 
        and summation either over the Log(OR) or OR.
        Returns naive polygenetic risk score (prs) as simple summation of odds ratios / log odds ratios (OR)
    """
    
    assert type(genes)==list
    assert type(csv)==str
    
    catalog=pd.read_csv(csv, delimiter = ',') #reading in the whole csv as a dataframe
    pat = '|'.join(r"{}".format(x) for x in genes) #compile search pattern for input genes of interest
    catalog_filtered = catalog[catalog["gene_name"].str.contains(pat,na=False)] #filtered catalog containing genes of interest only
    #print ("Catalog containing the inut genes only",catalog_filtered)
    #print ()
    pd_only_series=catalog_filtered["phewas phenotype"]=="Parkinson's disease" #filter catalog on pd phenotype only
    catalog_filtered_pd=catalog_filtered[pd_only_series]
    catalog_filtered_pd_pvalue_series=catalog_filtered_pd["p-value"]<p_threshold
    catalog_filtered_pd_pvalue=catalog_filtered_pd[catalog_filtered_pd_pvalue_series]
    
    #print (catalog_filtered_pd_pvalue)
    
    naive_prs=catalog_filtered_pd_pvalue["odds-ratio"].sum() #calculation of prs by choice of the user
    if log==True:
        return math.log10(naive_prs)
    else:
        return naive_prs


# In[2]:


# Mapping
rda_dic = pyreadr.read_r('LordickData/RiskSNPs.rda')
comb_df=rda_dic["comb"]
comb_df_PD=rda_dic["combPD"]


# In[4]:


comb_dic=comb_df.to_dict("list")
comb_df


# In[3]:


# first approach: Using the "neurommsig_pd_full_genes.csv" MechanismToGenes mapping to map the genes to SNPs in PheWAS
# PD associated MechanismToGenes dictionary 
genes=pd.read_csv('LordickData/neurommsig_pd_full_genes.csv')
MechanismToGene=genes.to_dict("list")


# In[16]:


print("Mechanisms Count:",len(MechanismToGene.keys()))
MechanismToGene


# In[57]:


count=0
genes_total=set()
for key in MechanismToGene.keys():
    values=[x for x in MechanismToGene[key] if type(x)!=float]
    genes_total.update(values)
    count+=len(values)
print ("Associated genes:",count)
print("total genes:",len(genes_total))


# In[58]:


# calculating the prs for each mechanism using the mapped genes only
count_noprs=0
count_prs=0
for mechanism in MechanismToGene.keys():
    print ("PD-associated mechanism: ", mechanism)
    prs_= prs(MechanismToGene[mechanism],"phewascatalog_full.csv",log=False)
    print("PRS: ",prs_)
    if prs_>0.0:
        count_prs+=1
    else:
        count_noprs+=1


# In[59]:


print ("Mechanisms with genes that can be associated to PD SNPs in PheWAS:", count_prs)
print ("Mechanisms with genes that can't be associated to PD SNPS:", count_noprs)


# Most gene lists can't be mapped to "Parkinson's Disease" SNPs in the PheWAS Database.
# So, the majority of all mechansims has a prs of zero.
# 
# # This is regardless of the eQTL,LD extented SNPs mapping
# 
# 
# Now lets check out the mapping:
# 

# In[21]:


# Second approach including the RiskSNP mapping from RiskSNPs.rda
rda_dic = pyreadr.read_r('LordickData/RiskSNPs.rda')
comb_df=rda_dic["comb"] # What's the difference between combPD and comb ? 
comb_PD_df=rda_dic["combPD"] 
print (comb_PD_df)


# # combPD data

# In[22]:


liste=comb_PD_df["rsID"].tolist()
liste_2=comb_PD_df["GENCODE_name"].tolist()
print ("SNPs in dataframe:",len(liste))
print ("SNPs in total:",len(set(liste)))
print ("Genes in dataframe:",len(liste_2))
print ("Genes in total:",len(set(liste_2)))


# In[23]:


#convert the mapping from the dataframe to a dictionary
GenetoSNP_PD={}
for x, y in zip(comb_PD_df['rsID'],comb_PD_df['GENCODE_name']):
    if y in GenetoSNP_PD:
        GenetoSNP_PD[y].append(x)
    else:
        GenetoSNP_PD[y]=[x]
print (GenetoSNP_PD)


# In[24]:


print(GenetoSNP_PD.keys())
print(len(GenetoSNP_PD.keys()))


# In[25]:


print (MechanismToGene)


# In[75]:


#pd.options.display.max_rows
#pd.set_option('display.max_rows', None)


# In[26]:


catalog=pd.read_csv("phewascatalog_full.csv",delimiter = ',') #reading in the whole csv as a dataframe


# In[27]:


#updating prs function to screen for snps in PheWAS

def prs_SNP_PD(gene,csv,log=False,p_threshold=0.05):
    
    #print ("SNPs:",GenetoSNP_PD[gene])
    #print("########################")
    if gene in GenetoSNP_PD:
        SNP_pat = '|'.join(r"{}".format(x) for x in GenetoSNP_PD[gene]) #compile search pattern for input genes of interest
        catalog_filtered = catalog[catalog["snp"].str.contains(SNP_pat,na=False)] #filtered catalog containing genes of interest only
        print ("Catalog containing the input snps only",catalog_filtered)
        print ("########################")

        catalog_filtered_pvalue_series=catalog_filtered["p-value"]<p_threshold
        catalog_filtered_pvalue=catalog_filtered[catalog_filtered_pvalue_series]

        #print (catalog_filtered_pd_pvalue)

        naive_prs=catalog_filtered_pvalue["odds-ratio"].sum() #calculation of prs by choice of the user
        
        if log==True:
            return math.log10(naive_prs)
        else:
            return naive_prs
    else:
        return False


# In[33]:


GenetoSNP_PD


# In[90]:


count=0
for gene in GenetoSNP_PD.keys():
    score=prs_SNP_PD(gene,catalog,log=False)
    if score>0.0:
        print (gene,":",score)
        count+=1
print("Total of genes with prs larger than 0 regardless of mechanism mapping:", count)


# Insgesamt 213 Gene, 61 mit PRS > 0 -> ~ 19 %

# Now: Mechanism -> Gene(s) -> SNPs

# In[28]:


#quickly remove nan values from gene lists to enhance efficiency
def remove_float_from_list(the_list):
    return [value for value in the_list if type(value)!=float]

MechanismToGene_noNans={x:remove_float_from_list(y) for x,y in MechanismToGene.items() }
print (MechanismToGene_noNans.keys())


# In[32]:


for i in MechanismToGene_noNans["Notch signaling subgraph"]:
    print (i)
    prs_SNP_PD(i,catalog)


# In[52]:


# Mechanism PRS Calculation
MechanismToScore=dict() #dictionary with gene to snps score as well as summated scores for each mechanism
for mechanism in MechanismToGene_noNans.keys():
    MechanismToScore[mechanism]=[[],[0]] 
    for gene in MechanismToGene_noNans[mechanism]:
        gene,score=gene,prs_SNP_PD(gene,catalog,log=False)
        MechanismToScore[mechanism][0].append((gene,score))
        if score!=False:
            MechanismToScore[mechanism][1]+=score


# In[53]:


MechanismToScore


# Almost most of the genes appear to have no prs because it's not in the mapping!

# In[81]:


# actual fraction of PD related Genes of the Mechanism-gene mapping that are actually in the Gene-SNP mapping
gene_list_ofallmechanisms=[y for x in list(MechanismToGene_noNans.values()) for y in x]
gene_set_ofallmechanisms=list(set(gene_list_ofallmechanisms))
count=0
for i in gene_set_ofallmechanisms:
    if i in GenetoSNP_PD:
        count+=1
print("Number of PD related genes:",len(gene_set_ofallmechanisms))
print("Number of PD related genes occuring in the GeneToSNP mapping:", count)
print ("Fraction:", count/len(gene_set_ofallmechanisms))


# Problem: Only ~ 5 % of the PD related genes from the NeuroMMSIG thing occure in the GeneToSNP mapping

# In[57]:


#more summary statistics: How many Mechanisms actually have useful scores?
count=0
for mechanism in MechanismToScore:
    if MechanismToScore[mechanism][1][0]>0:
        count+=1
print("Fraction of PD Mechanism with a PRS > 0:",count/len(MechanismToScore))


# Big overlap! very few genes map to almost half of the PD Mechanisms! 

# # comb data

# In[87]:


# comb data obviously lesser dataframe size but more snps and genes in total
liste=comb_df["rsID"].tolist()
liste_2=comb_df["GENCODE_name"].tolist()
print ("SNPs in dataframe:",len(liste))
print ("SNPs in total:",len(set(liste)))
print ("Genes in dataframe:",len(liste_2))
print ("Genes in total:",len(set(liste_2)))


# In[19]:


#updating prs_SNO function for the comb genes-snps
def prs_SNP(gene,csv,log=False,p_threshold=0.05):
    
    #print ("SNPs:",GenetoSNP_PD[gene])
    #print("########################")
    if gene in GenetoSNP:
        #catalog=pd.read_csv(csv, delimiter = ',') #reading in the whole csv as a dataframe
        SNP_pat = '|'.join(r"{}".format(x) for x in GenetoSNP[gene]) #compile search pattern for input genes of interest

        catalog_filtered = catalog[catalog["snp"].str.contains(SNP_pat,na=False)] #filtered catalog containing genes of interest only
        #print ("Catalog containing the input snps only",catalog_filtered)
        #print ("########################")
        #pd_only_series=catalog_filtered["phewas phenotype"]=="Parkinson's disease" #filter catalog on pd phenotype only
        #catalog_filtered_pd=catalog_filtered[pd_only_series]
        catalog_filtered_pvalue_series=catalog_filtered["p-value"]<p_threshold
        catalog_filtered_pvalue=catalog_filtered[catalog_filtered_pvalue_series]

        #print (catalog_filtered_pd_pvalue)

        naive_prs=catalog_filtered_pvalue["odds-ratio"].sum() #calculation of prs by choice of the user
        if log==True:
            return math.log10(naive_prs)
        else:
            return naive_prs
    else:
        return False


# In[95]:


#convert the mapping from the dataframe to a dictionary
GenetoSNP={}
for x, y in zip(comb_df['rsID'],comb_df['GENCODE_name']):
    if y in GenetoSNP:
        GenetoSNP[y].append(x)
    else:
        GenetoSNP[y]=[x]
print (GenetoSNP)


# In[91]:


count=0
for gene in GenetoSNP.keys():
    score=prs_SNP(gene,catalog,log=False)
    if score>0.0:
        print (gene,":",score)
        count+=1
print("Total of genes with prs larger than 0 regardless of mechanism mapping:", count)


# In[92]:


print ("Fraction of genes with prs larger than 0 regardless of mechansim mapping:",21/616)


# Now: Mechanism -> Gene(s) -> SNPs

# In[98]:


# Mechanism PRS Calculation
for mechanism in MechanismToGene_noNans.keys():
    MechanismToScore[mechanism]=[[],[0]] 
    for gene in MechanismToGene_noNans[mechanism]:
        gene,score=gene,prs_SNP(gene,catalog,log=False)
        MechanismToScore[mechanism][0].append((gene,score))
        if score!=False:
            MechanismToScore[mechanism][1]+=score


# In[99]:


MechanismToScore


# In[100]:


count=0
for i in gene_set_ofallmechanisms:
    if i in GenetoSNP:
        count+=1
print("Number of PD related genes:",len(gene_set_ofallmechanisms))
print("Number of PD related genes occuring in the GeneToSNP mapping:", count)
print ("Fraction:", count/len(gene_set_ofallmechanisms))


# In[101]:


#more summary statistics: How many Mechanisms actually have useful scores?
count=0
for mechanism in MechanismToScore:
    if MechanismToScore[mechanism][1][0]>0:
        count+=1
print("Fraction of PD Mechanism with a PRS > 0:",count/len(MechanismToScore))


# In[35]:


####################################################################################################


# In[104]:


####### To do: compare original SNPs (~ 3700) to mapping SNPS maybe #######


# In[54]:


snps_3700=pd.read_csv('LordickData/PPMI_combPD_pass_CADD.csv')
snpsIDs_3700=snps_3700["ID"].tolist()
print ("count:", len(snps_3700["ID"]))


# In[64]:


mapping_snps=list(GenetoSNP_PD.values())
mapping_snps = [item for items in mapping_snps for item in items]
mapping_snps_set=set(mapping_snps)


# In[66]:


print("count SNPS from Mapping file:",len(mapping_snps_set))


# In[58]:


# -> more SNPs in mapping than in SNP origin database?!


# In[68]:


# overlap 
len(mapping_snps_set.intersection(snpsIDs_3700))


# All SNPs from the Origin (3742) fall into the mapping SNPS (RiskSnps). Hence: extra 500 SNPs are in the RiskSNP dataset for combPD

# In[75]:


#check quickly how many snps we have have are actually in the pheWAS catalog
def check4SNPs(snpslist,catalog):
    
    #catalog=pd.read_csv(csv, delimiter = ',') #reading in the whole csv as a dataframe
    SNP_pat = '|'.join(r"{}".format(x) for x in snpslist) #compile search pattern for input genes of interest
    catalog_filtered = catalog[catalog["snp"].str.contains(SNP_pat,na=False)] #filtered catalog containing genes of interest only
    print ("count of input SNPs:",len(snpslist))
    print ("count of input SNPs found in PheWAS:",len(catalog_filtered))
    
    return ("SNPs actually missing in PheWAS: "+str(len(snpslist)-len(catalog_filtered)))


# In[74]:


#check for origin SNPs 3700
check4SNPs(snpsIDs_3700,catalog)


# In[76]:


#check for mapping snps
check4SNPs(list(mapping_snps_set),catalog)


# In[ ]:


# -> 54 of the additional 500 SNPs in the mapping are also found in the db. No big difference at all. 

