import twobitreader
import numpy as np
import pandas as pd
import os

genome = twobitreader.TwoBitFile("hg19.2bit")

def getseq_from_genome(from_genome,chromosomenum, start_pos, end_pos):
    #from_genome = twobitreader.TwoBitFile("hg19.2bit")
    chromosome = from_genome[chromosomenum]
    seq = chromosome.get_slice(start_pos-1, end_pos)
    return seq.upper()

vcf_filename = "1000G.EUR.QC.plink.simple.vcf"
usecolumns_list=[0,1,2,3,4]
column_renames_list=['CHR','POS','ID','REF','ALT']
outfolder="genome_database"
window_size=30
os.mkdir(outfolder)
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
    



