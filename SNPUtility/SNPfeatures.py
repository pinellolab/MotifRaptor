import numpy as np
import pandas as pd
import os
import scipy.stats as stats
import pybedtools
import re
import multiprocessing
from functools import partial
import pyBigWig

def SNPExtractFeatureFromBigwig(featurefile,SNPbedfile):
    SNP_df=pd.read_csv(SNPbedfile, sep="\t", header=None)
    bw = pyBigWig.open(featurefile)
    valuelist=[]
    for index in SNP_df.index:
        cols0=SNP_df.iloc[index][0]
        cols1=SNP_df.iloc[index][1]
        cols2=SNP_df.iloc[index][2]
        vals = bw.values(cols0,cols1,cols2)[0]
        if np.isnan(vals):
            valuelist.append(0)
        else:
            valuelist.append(vals)
    bw.close()
    n=len(SNP_df.columns)
    SNP_df[n]=valuelist
    return SNP_df


def SNPExtractFeatureFromBedfile(featurefile,SNPbedfile):
    a = pybedtools.BedTool(SNPbedfile)
    a=a.sort()
    b=pybedtools.BedTool(featurefile)
    #b=b.sort()
    c = a.closest(b, d=True)
    c.moveto("tempfeatures.bed")
    SNP_df=pd.read_csv("tempfeatures.bed", sep="\t", header=None)
    SNP_df_sub=SNP_df.loc[SNP_df[8]==0]
    SNP_df_2=SNP_df_sub[[0,1,2,3,7]].copy()
    return SNP_df_2


def SNPfeatures_main(featurefolder, SNPbedfile, SNPfeaturedf_outfile, SNPnormalize_list):
    filenames=os.listdir(featurefolder)
    feature_df_master=pd.read_csv(SNPbedfile, sep="\t", header=None)
    feature_df_master.columns=['chr','start','end','snpid']
    for filename in filenames:
    #filename='CADD.bw'
    #if filename=='CADD.bw':
        print("Processing..."+filename+"!!!")
        extname=filename.split(".")[-1]
        fileprefix=filename.split(".")[0]
        if extname=="bw":
            featurefile=os.path.join(featurefolder,filename)
            feature_df=SNPExtractFeatureFromBigwig(featurefile,SNPbedfile)
            feature_df.columns=['chr','start','end','snpid',fileprefix]
            feature_df_master = pd.merge(feature_df_master,
                     feature_df[['snpid',fileprefix]],
                     on='snpid',
                     how='left')
        elif extname=="bed":
            featurefile=os.path.join(featurefolder,filename)
            feature_df=SNPExtractFeatureFromBedfile(featurefile,SNPbedfile)
            feature_df.columns=['chr','start','end','snpid',fileprefix]
            feature_df_master = pd.merge(feature_df_master,
                     feature_df[['snpid',fileprefix]],
                     on='snpid',
                     how='left')

    feature_df_master=feature_df_master.fillna(0)
    feature_df_master_sub=feature_df_master[feature_df_master['snpid'].isin(SNPnormalize_list)]
    column_max=[]
    column_min=[]
    for column_num in range(4,feature_df_master.shape[1]):
        column_name=feature_df_master.columns[column_num]
        a=np.array(feature_df_master_sub[column_name])
        column_max.append(np.percentile(a,99))
        column_min.append(np.percentile(a,1))

    maxrow=['0','0','0','max']+column_max
    minrow=['0','0','0','min']+column_min
    df2 = pd.DataFrame([maxrow,minrow],columns=feature_df_master.columns)
    feature_df_master=pd.concat([df2,feature_df_master])
    feature_df_master.to_csv(SNPfeaturedf_outfile,index=None, header=True, sep="\t")
