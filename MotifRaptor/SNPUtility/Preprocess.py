import twobitreader
import numpy as np
import pandas as pd
import os
import scipy.stats as st

def getseq_from_genome(from_genome,chromosomenum, start_pos, end_pos):
    #from_genome = twobitreader.TwoBitFile("hg19.2bit")
    chromosome = from_genome[chromosomenum]
    seq = chromosome.get_slice(start_pos-1, end_pos)
    return seq.upper()

def preprocess_vcf_for_scan(genome, vcf_filename, usecolumns_list, column_renames_list, outfolder, window_size): 
    #genome = twobitreader.TwoBitFile("hg19.2bit")
    genome = twobitreader.TwoBitFile(genome)
    #vcf_filename = "1000G.EUR.QC.plink.simple.vcf"
    #usecolumns_list=[0,1,2,3,4]
    #column_renames_list=['CHR','POS','ID','REF','ALT']
    #outfolder="genome_database"
    #window_size=30
    #os.mkdir(outfolder)
    df_snps_seq=pd.read_table(vcf_filename,sep='\t',header=0,
                      usecols=usecolumns_list, 
                      names=column_renames_list, 
                      low_memory=False, comment='#')
    df_snps_seq=df_snps_seq.set_index('ID')

    df_snps_seq['UID']=df_snps_seq.index+":"+df_snps_seq['REF']+":"+df_snps_seq['ALT']

    seqfiledict=dict()
    posfiledict=dict()

    chromosomeset=set(df_snps_seq['CHR'])

    for chromosome in chromosomeset:
        chromosomenum="chr"+str(chromosome).strip()
        df_snp_out_sub=df_snps_seq[df_snps_seq['CHR']==chromosome]['UID']
        df_snp_out_sub.to_csv(os.path.join(outfolder, chromosomenum+".uid"),sep='\t',index=True,header=True)

    for chromosome in chromosomeset:
        chromosomenum="chr"+str(chromosome).strip()
        seqfilename=os.path.join(outfolder,chromosomenum+"_SEQ_REF_ALT.txt")
        posfilename=os.path.join(outfolder,chromosomenum+"_POS_REF_ALT.txt")
        seqfile=open(seqfilename,"w")
        posfile=open(posfilename,"w")
        seqfiledict[seqfilename]=seqfile
        posfiledict[posfilename]=posfile
        #posfile.write("START_POS"+"\t"+"END_POS"+"\t"+"SNP_SEQ"+"\t"+"ref_or_alt"+"\n")

    snpcount=0
    linecount=0
    basepaircount=0
    previouschromosome="chr0"
    for snpid in df_snps_seq.index:
        currentsnp=df_snps_seq.loc[snpid]
        chromosomenum="chr"+str(currentsnp['CHR']).strip()
        if(chromosomenum!=previouschromosome):
            print("Current number of SNP :"+str(snpcount))
            print("Current line number: "+str(linecount))
            print("Current base pair number: "+str(basepaircount))
            previouschromosome=chromosomenum
        seqfilename=os.path.join(outfolder,chromosomenum+"_SEQ_REF_ALT.txt")
        posfilename=os.path.join(outfolder,chromosomenum+"_POS_REF_ALT.txt")
        seqfile=seqfiledict[seqfilename]
        posfile=posfiledict[posfilename]
        snppos=int(currentsnp['POS'])
        start_pos=int(currentsnp['POS'])-30
        end_pos=int(currentsnp['POS'])+30
        snpchar=getseq_from_genome(genome,chromosomenum, snppos, snppos)
        seq=getseq_from_genome(genome,chromosomenum, start_pos, end_pos)
        ref=str(currentsnp['REF']).strip().upper()
        alt=str(currentsnp['ALT']).strip().upper()
        ref_seq=(seq[0:30]+ref+seq[31:61]).upper()
        alt_seq=(seq[0:30]+alt+seq[31:61]).upper()
        if(len(ref_seq)!=61 or len(alt_seq)!=61) :
            print("something_wrong!!!")
            print("ref_seq: "+ref_seq)
            print("alt_seq: "+alt_seq)
            print("ID: "+snpid)
            print("SNP :"+snpchar)
            break
        #if 'N' in ref_seq or 'N' in alt_seq:
        #    continue
        seqfile.write(ref_seq)
        seqfile.write(alt_seq)
        basepaircount=basepaircount+len(ref_seq)+len(alt_seq)
        #refposstring=str(start_pos)+"\t"+str(end_pos)+"\t"+ref+"\t"+"ref"+"\n"
        #altposstring=str(start_pos)+"\t"+str(end_pos)+"\t"+alt+"\t"+"alt"+"\n"
        refposstring=str(start_pos)+"\n"
        altposstring=str(start_pos)+"\n"
        posfile.write(refposstring)
        posfile.write(altposstring)
        linecount=linecount+2
        snpcount=snpcount+1
    print("Total number of SNP :"+str(snpcount))
    print("Total line number: "+str(linecount))

    for f in seqfiledict.keys():
        seqfiledict[f].close()
    for f in posfiledict.keys():
        posfiledict[f].close()
    
def preprocess_summary_statistics(sumstatsfile):
    #sumstatsfile="RA_GWASmeta_TransEthnic_v2.txt"
    SNP_sumstat_dataframe=pd.read_csv(sumstatsfile,sep='\t')
    SNP_sumstat_dataframe=SNP_sumstat_dataframe[SNP_sumstat_dataframe.columns[[0,1,2,3,4,6]]]
    SNP_sumstat_dataframe.columns=['ID','CHR','POS','A1','A2','P-val']
    p=5E-8
    hit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']<=p]
    nonhit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']>p]
    hit_SNP_df_sub=hit_SNP_df[['ID','CHR','POS']]
    hit_SNP_df_sub.to_csv("hitSNP_list.txt", sep='\t',index=None, header=True)
    nonhit_SNP_df_sub=nonhit_SNP_df[['ID','CHR','POS']]
    nonhit_SNP_df_sub.to_csv("nonhitSNP_list.txt", sep='\t',index=None, header=True)
    hit_SNP_cvf_sub=hit_SNP_df[['ID','CHR','POS','A1','A2']]
    hit_SNP_vcf_sub.to_csv("hitSNP_list.vcf", sep='\t',index=None, header=True)

def preprocess_summary_statistics_general(sumstatsfile,score_type,pthreshold,column_num_list):
    #sumstatsfile="RA_GWASmeta_TransEthnic_v2.txt"
    SNP_sumstat_dataframe=pd.read_csv(sumstatsfile,sep='\t')
    SNP_sumstat_dataframe=SNP_sumstat_dataframe[SNP_sumstat_dataframe.columns[column_num_list]]
    if score_type=="pvalue":
        SNP_sumstat_dataframe.columns=['ID','CHR','POS','A1','A2','P-val']
        #p=5E-8
        p = pthreshold
        hit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']<=p]
        nonhit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']>p]
        hit_SNP_df_sub=hit_SNP_df[['ID','CHR','POS']]
        hit_SNP_df_sub.to_csv("hitSNP_list.txt", sep='\t',index=None, header=True)
        nonhit_SNP_df_sub=nonhit_SNP_df[['ID','CHR','POS']]
        nonhit_SNP_df_sub.to_csv("nonhitSNP_list.txt", sep='\t',index=None, header=True)
        hit_SNP_vcf_sub=hit_SNP_df[['ID','CHR','POS','A1','A2']]
        hit_SNP_vcf_sub.to_csv("hitSNP_list.vcf", sep='\t',index=None, header=True)
    elif score_type=="zscore":
        zcutoff=st.norm.ppf(pthreshold)
        SNP_sumstat_dataframe.columns=['ID','CHR','POS','A1','A2','Z']
        hit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['Z']>np.abs(zcutoff)]
        nonhit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['Z']<=np.abs(zcutoff)]

        hit_SNP_df_sub=hit_SNP_df[['ID','CHR','POS']]
        hit_SNP_df_sub.to_csv("hitSNP_list.txt", sep='\t',index=None, header=True)
        nonhit_SNP_df_sub=nonhit_SNP_df[['ID','CHR','POS']]
        nonhit_SNP_df_sub.to_csv("nonhitSNP_list.txt", sep='\t',index=None, header=True)
        hit_SNP_vcf_sub=hit_SNP_df[['ID','CHR','POS','A1','A2']]
        hit_SNP_vcf_sub.to_csv("hitSNP_list.vcf", sep='\t',index=None, header=True)
    elif score_type=="chisquare":
        SNP_sumstat_dataframe.columns=['ID','CHR','POS','A1','A2','CHI']
        p_list=[]
        for ind in SNP_sumstat_dataframe.index:
            chiscore = SNP_sumstat_dataframe.loc[ind]['CHI']
            pvalue=1-st.chi2.cdf(float(chiscore), 1)
            p_list.append(pvalue)
        SNP_sumstat_dataframe['P-val']=p_list
        p = pthreshold
        hit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']<=p]
        nonhit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']>p]
        hit_SNP_df_sub=hit_SNP_df[['ID','CHR','POS']]
        hit_SNP_df_sub.to_csv("hitSNP_list.txt", sep='\t',index=None, header=True)
        nonhit_SNP_df_sub=nonhit_SNP_df[['ID','CHR','POS']]
        nonhit_SNP_df_sub.to_csv("nonhitSNP_list.txt", sep='\t',index=None, header=True)
        hit_SNP_vcf_sub=hit_SNP_df[['ID','CHR','POS','A1','A2']]
        hit_SNP_vcf_sub.to_csv("hitSNP_list.vcf", sep='\t',index=None, header=True)


def preprocess_from_ukbb_v3(ukbb_tsv_filename,pthreshold):
    #sumstatsfile="30780_raw.gwas.imputed_v3.both_sexes.tsv"
    sumstatsfile = ukbb_tsv_filename
    SNP_sumstat_dataframe=pd.read_csv(sumstatsfile,sep='\t')

    SNPID_list=[]
    CHR_list=[]
    POS_list=[]
    A1_list=[]
    A2_list=[]
    pvalue_list=[]

    for index in SNP_sumstat_dataframe.index:
        thissnp=SNP_sumstat_dataframe.loc[index]
        variant_element=thissnp['variant'].split(":")
        #chromnum="chr"+str(variant_element[0])
        chromnum=str(variant_element[0])
        if chromnum=='X' or chromnum=='Y' or chromnum=='M':
            continue
        position_based_one=str(variant_element[1])
        minorallele=thissnp['minor_allele']
        major_allele=''
        minor_allele=''
        if minorallele==str(variant_element[2]):
            major_allele=str(variant_element[3])
            minor_allele=str(variant_element[2])
        elif minorallele==str(variant_element[3]):
            major_allele=str(variant_element[2])
            minor_allele=str(variant_element[3])
        snpid=chromnum+"_"+position_based_one+"_"+major_allele+"_"+minor_allele
        SNPID_list.append(snpid)
        CHR_list.append(chromnum)
        POS_list.append(position_based_one)
        A1_list.append(major_allele)
        A2_list.append(minor_allele)
        pvalue_list.append(float(thissnp['pval']))

    outputdf=pd.DataFrame(columns=['ID','CHR','POS','A1','A2','P-val'])
    outputdf['ID']=SNPID_list
    outputdf['CHR']=CHR_list
    outputdf['POS']=POS_list
    outputdf['A1']=A1_list
    outputdf['A2']=A2_list
    outputdf['P-val']=pvalue_list

    SNP_sumstat_dataframe=outputdf
    #p=5E-8
    p = pthreshold
    hit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']<=p]
    nonhit_SNP_df=SNP_sumstat_dataframe[SNP_sumstat_dataframe['P-val']>p]
    hit_SNP_df_sub=hit_SNP_df[['ID','CHR','POS']]
    hit_SNP_df_sub.to_csv("hitSNP_list.txt", sep='\t',index=None, header=True)
    nonhit_SNP_df_sub=nonhit_SNP_df[['ID','CHR','POS']]
    nonhit_SNP_df_sub.to_csv("nonhitSNP_list.txt", sep='\t',index=None, header=True)
    hit_SNP_vcf_sub=hit_SNP_df[['ID','CHR','POS','A1','A2']]
    hit_SNP_vcf_sub.to_csv("hitSNP_list.vcf", sep='\t',index=None, header=True)

