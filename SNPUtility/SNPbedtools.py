import numpy as np
import pandas as pd
import os
import scipy.stats as stats
import pybedtools
import re
import multiprocessing
from functools import partial

import matplotlib.pyplot as plt

def generateSNPBed(SNPfilename,SNPbedfilename):
    output_hitdataframe=pd.read_csv(SNPfilename,sep='\t')
    SNP_score_dataframe_sub=output_hitdataframe
    SNP_score_dataframe_sub=SNP_score_dataframe_sub.set_index('ID')
    output_bed=SNPbedfilename+"_tempSNP_unsorted.bed"
    output_sorted_bed=SNPbedfilename

    outputfile=open(output_bed, "w")
    for snpid in SNP_score_dataframe_sub.index:
        thissnp=SNP_score_dataframe_sub.loc[snpid]
        outputstring="chr"+str(thissnp.CHR)+"\t"+str(thissnp.POS-1)+"\t"+str(thissnp.POS)+"\t"+snpid
        outputfile.write(outputstring+"\n")
    outputfile.close()
    sortcommand="sort -k1,1 -k2,2n "+output_bed+" > "+output_sorted_bed
    os.system(sortcommand)

def generateSNPBed_fast(SNPfilename,SNPbedfilename):
    output_hitdataframe=pd.read_csv(SNPfilename,sep='\t',header=0)
    SNP_score_dataframe_sub=output_hitdataframe
    SNP_score_dataframe_sub=output_hitdataframe
    SNP_score_dataframe_sub['POS1']=SNP_score_dataframe_sub['POS']-1
    chromlist=SNP_score_dataframe_sub['CHR']
    SNP_score_dataframe_sub['CHR1']=["chr" + str(c) for c in chromlist]
    output_bed=SNPbedfilename+"_tempSNP_unsorted.bed"
    

    SNP_score_dataframe_sub_out=SNP_score_dataframe_sub[['CHR1','POS1','POS','ID']]
    SNP_score_dataframe_sub_out.to_csv(output_bed,sep="\t",header=None, index=None)
    
    output_sorted_bed=SNPbedfilename
    sortcommand="sort -k1,1 -k2,2n "+output_bed+" > "+output_sorted_bed
    os.system(sortcommand)


def overlapSNPandBed(SNPsortedbed, Bedfile, outputBedfile):
    output_sorted_bed=SNPsortedbed
    a_union_bg = pybedtools.BedTool(Bedfile)#DHS_process/union_DHS.bed"
    b = pybedtools.BedTool(output_sorted_bed)
    (b+a_union_bg).moveto(outputBedfile)#"hitSNP_DHSunion_list.txt"


