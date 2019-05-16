import twobitreader
import numpy as np
import pandas as pd
import os
import sys
import argparse
import multiprocessing as mp
from functools import partial

def getseq_from_genome(from_genome, chromosomenum, start_pos, end_pos):    
    chromosome = from_genome[chromosomenum]
    seq = chromosome.get_slice(start_pos-1, end_pos)
    return seq.upper()

def getseq_from_genome_2(from_genome, chromosomenum, start_pos, end_pos, mask_repetitive=False):   
    if mask_repetitive:
        seq= ''.join([mask(c) for c in from_genome[chromosomenum][start_pos-1:end_pos]])
    else:
        seq= from_genome[chromosomenum][start_pos-1:end_pos]
    return seq.upper()



def extract_sequence_ref_alt(genome,window_size,snp):
    #print(str(snp.ID)+" "+str(snp.REF)+" "+str(snp.ALT))
    ref_snp=snp.REF.replace('-','').strip()
    alt_snp=snp.ALT.replace('-','').strip()
    ref_snp_len=len(ref_snp)
    alt_snp_len=len(alt_snp)
    pos_snp=int(snp.POS)
    start_pos=int(pos_snp-window_size)
    end_pos=start_pos+int(2*window_size+ref_snp_len)
    chromnum=str(int(snp.CHR))
    chromosomenum="chr"+str(chromnum).strip()
    seq=getseq_from_genome(genome,chromosomenum, start_pos, end_pos)
    pos= int(pos_snp-start_pos)

    seq_ref= ''.join([seq[:pos],ref_snp,seq[pos+ref_snp_len:]])
    seq_alt= ''.join([seq[:pos],alt_snp,seq[pos+ref_snp_len:]])
    rawrefseq=seq[pos:pos+ref_snp_len].upper()
    #if rawrefseq!=ref_snp:
    #    print("WRONG SNP: "+snp.ID)
    snp['seq_ref']=seq_ref
    snp['seq_alt']=seq_alt
    snp['pos_start_ref']=pos
    snp['pos_end_ref']=pos+ref_snp_len-1
    snp['pos_start_alt']=pos
    snp['pos_end_alt']=pos+alt_snp_len-1
    
    return snp

def process_chunk(genome, window_size, df):
    #genome = twobitreader.TwoBitFile("hg19.2bit")
    compute_ref_alt=partial(extract_sequence_ref_alt,genome,window_size)
    return df.apply(compute_ref_alt, axis=1)

print("python VCFtoSeq.py vcf_filename genome_filename usecolumns column_renames onesidewindowsize numberofthreads outputpicklefile")
print("Example python VCFToSeq.py 1000G.EUR.QC.plink.simple.vcf hg19.2bit 0,1,2,3,4 CHR,POS,ID,REF,ALT 30 4 hitpickle")
print("genome_filename: must be a 2bit file")
print("usecolumns: start from 0")
print("column_renames: rename the columns, aligned with usecolumns, must at least contains ID,CHR,POS,REF,ALT")
print("hitpicke: output will be hitpickle.pickle and hitpickle.df.txt")

vcf_filename=sys.argv[1]
genome_filename=sys.argv[2]
usecolumnsstring=sys.argv[3].strip().split(',')
usecolumns_list=[int(ai) for ai in usecolumnsstring]

column_renames_list=sys.argv[4].strip().split(',')
window_size=int(sys.argv[5])
numberofthreads=int(sys.argv[6])
outputpicklefile=sys.argv[7]

df_snps=pd.read_table(vcf_filename,sep='\t',header=0,
                      usecols=usecolumns_list, 
                      names=column_renames_list, 
                      low_memory=False, comment='#')
#df_snps.set_index('ID', inplace=True)
#df_snps.to_csv(outputpicklefile+".df.txt",sep='\t')
#window_size=30
mask_repetitive=False
genome = twobitreader.TwoBitFile(genome_filename)
process_chunk_partial=partial(process_chunk, genome, window_size)
p = mp.Pool(processes=numberofthreads)
pool_results = p.map(process_chunk_partial,np.array_split(df_snps,numberofthreads))
p.close()
p.join()

# merging parts processed by different processes
df_snps_seq = pd.concat(pool_results, axis=0)

df_snps_seq.set_index('ID', inplace=True)
df_snps_seq.to_pickle(outputpicklefile+".pickle")
df_snps_seq.to_csv(outputpicklefile+".df.txt",sep='\t',index=True)
