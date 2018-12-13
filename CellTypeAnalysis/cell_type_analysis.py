import numpy as np
import pandas as pd
import os
import scipy.stats as stats
import pybedtools
import re
import multiprocessing
from functools import partial
from statsmodels.sandbox.stats.multicomp import multipletests

import seaborn as sns
sns.set_style('white')
sns.set_context('poster')

def celltype_analysis(target_SNP_df, bg_SNP_df, output_sorted_BedTool_hit,output_sorted_BedTool_nonhit, output_dir_motif, DHSdir, celltypetable,celltype):        
    fileaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='DNase-seq')]['File accession'].iloc[0]
    targetbedfile=os.path.join(DHSdir, fileaccession+'_target.bed')
    bgbedfile=os.path.join(DHSdir, fileaccession+'_balanced_bg.bed')
    a = pybedtools.BedTool(targetbedfile)
    a_bg = pybedtools.BedTool(bgbedfile)
    targetresultbedfile=os.path.join(output_dir_motif,celltype+'.result.bed')
    bgresultbedfile=os.path.join(output_dir_motif,celltype+'.result_bg.bed')
    (output_sorted_BedTool_hit+a).moveto(targetresultbedfile)
    #(output_sorted_BedTool_nonhit+a_bg).moveto(bgresultbedfile)
    (output_sorted_BedTool_nonhit+a).moveto(bgresultbedfile)
    SNP1=pd.read_csv(targetresultbedfile,sep='\t',header=None)
    SNP1.columns=['CHR','START','STOP','ID']
    SNP1=SNP1.set_index('ID')
    
    SNP_hitscore_fg_filtered=SNP1
    SNP2=pd.read_csv(bgresultbedfile,sep='\t',header=None)
    SNP2.columns=['CHR','START','STOP','ID']
    SNP2=SNP2.set_index('ID')
    
    SNP_nonhitscore_fg_filtered=SNP2
    hit_target_index=set(SNP_hitscore_fg_filtered.index)
    nonhit_target_index=set(SNP_nonhitscore_fg_filtered.index)
    q=len(hit_target_index) #1
    #m=len(set(target_SNP_df.index)) #3
    m=output_sorted_BedTool_hit.count()
    p=len(nonhit_target_index) #2
    #n=len(set(bg_SNP_df.index)) #4
    n=output_sorted_BedTool_nonhit.count()
    oddsratio, pvalue = stats.fisher_exact([[q, m-q], [p, n-p]])
    #oddsratio, pvalue = stats.fisher_exact([[q, p], [m, n]])
    outputfilename=os.path.join(output_dir_motif, celltype+".pvalue")
    outputfile=open(outputfilename,"w")
    outputstring=celltype+"\t"+fileaccession+"\t"+str(q)+"\t"+str(m)+"\t"+str(p)+"\t"+str(n)+"\t"+str(oddsratio)+"\t"+str(pvalue)
    outputfile.write(outputstring+"\n")
    outputfile.close()
    print(celltype+" Done!")

def run_celltype_analysis(target_SNP_df, bg_SNP_df, output_sorted_bed_hits,output_sorted_bed_nonhits, output_dir_motif, DHSdir, celltypetable,celltypelist, numberofthreads):
    b_hits = pybedtools.BedTool(output_sorted_bed_hits)
    b_nonhits = pybedtools.BedTool(output_sorted_bed_nonhits)
    p=multiprocessing.Pool(processes=numberofthreads)
    whole_analysis_func=partial(celltype_analysis, target_SNP_df, bg_SNP_df, b_hits, b_nonhits,output_dir_motif, DHSdir, celltypetable)
    p.map(whole_analysis_func, celltypelist)
    p.close()    

def run_plot_figure(pdffilename,result_pvalue):
    result_pvalue=result_pvalue.sort_values(by=[7],ascending=False)
    plt.rc('pdf', fonttype=42)
    plt.figure(figsize=(10,20))
    ybars=result_pvalue[0]
    y_pos = np.arange(len(ybars))
    xbars=-np.log10(result_pvalue[8].tolist())
    labels=["{0:.2f}".format(round(float(f),2)) for f in result_pvalue[6]]
    plt.barh(y_pos, xbars,color='blue')
    plt.plot([-np.log10(0.05),-np.log10(0.05)],[-100,100],'r--')
    plt.ylim([-2,84])
    plt.yticks(y_pos, ybars)
    plt.xlabel('-log10 corrected p value')

    for i in range(len(y_pos)):
        plt.text(x = xbars[i]+0.1 , y = y_pos[i]-0.3, s = labels[i], size = 9)
    
    plt.savefig(pdffilename, bbox_inches='tight', pad_inches=0)



def main():
    output_beds_dir="./testcelltype"
    celltypelistfile="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_08_new_data_structure/celltypelist.txt"
    celltypemappingfile="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_07_ENCODE_ROADMAP/download_meta_tier123.txt"
    #DHSdir="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_07_ENCODE_ROADMAP/DHS/"
    DHSdir="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_07_ENCODE_ROADMAP/DHS_less_20/"
    numberofthreads=4
    p=5E-8

    motif_id="bedfiles"
    output_dir_motif=os.path.join(output_beds_dir, motif_id)
    if not os.path.exists(output_dir_motif):
        os.mkdir(output_dir_motif)




    output_sorted_bed_hits=os.path.join(output_dir,"hitSNP_DHSunion_list.bed")
    output_sorted_bed_nonhits=os.path.join(output_dir,"nonhitSNP_DHSunion_list.bed")

    celltypelist=list(pd.read_csv(celltypelistfile,sep='\t',header=None)[0])
    celltypetable=pd.read_csv(celltypemappingfile,sep='\t')
    run_celltype_analysis(target_SNP_df, bg_SNP_df, output_sorted_bed_hits,output_sorted_bed_nonhits, output_dir_motif, DHSdir, celltypetable,celltypelist, numberofthreads)


    #!cat ./testcelltype/bedfiles/*.pvalue > ./testcelltype/all.pvalue
    from statsmodels.sandbox.stats.multicomp import multipletests
    result_pvalue=pd.read_csv("testcelltype/all.pvalue",sep="\t",header=None)
    result_pvalue=result_pvalue.sort_values(by=[6],ascending=False)
    pvalues=result_pvalue[7]
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]
    result_pvalue[8]=qvalue
    result_pvalue=result_pvalue.sort_values(by=[7],ascending=True)

    run_plot_figure('cell_type_plot.pdf', result_pvalue)
