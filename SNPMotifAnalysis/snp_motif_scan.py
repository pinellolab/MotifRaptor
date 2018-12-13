import twobitreader
import numpy as np
import pandas as pd
import os
import argparse
import multiprocessing
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

def generate_pfm_matrix(motif_pfm_filename, pseudo_count, background=None):
    motif_pfm_matrix=pd.read_csv(motif_pfm_filename, header=None,delim_whitespace=True)
    motif_pfm_matrix.index=['A','C','G','T']
    motif_pfm_matrix=motif_pfm_matrix.astype('float64')
    if background==None:
        bg = {'A':0.2952,'C':0.2048,'G':0.2048,'T':0.2952}
    else:
        bg=background
    for key in bg.keys():
        bg[key]=np.log(bg[key])
    motif_pfm_matrix=motif_pfm_matrix+pseudo_count
    pfm_sum=np.sum(motif_pfm_matrix,axis=0)
    motif_pfm_matrix_2=np.log(motif_pfm_matrix/pfm_sum)
    for nucleotide in motif_pfm_matrix_2.index:
        for column in motif_pfm_matrix_2.columns:
            motif_pfm_matrix_2.loc[nucleotide,column]=motif_pfm_matrix_2.loc[nucleotide,column]-bg[nucleotide]
    return motif_pfm_matrix_2

def sequence_motif_score(sequence, motif_pfm_matrix):
    score=0
    column=0
    for nucleotide in sequence:
        score = score + motif_pfm_matrix.loc[nucleotide,column]
        column=column+1
    return score

def generate_disrupt_dict(motif_pfm_matrix):
    disrupt_dict=dict()
    #max_disrupt_score_list=[]
    for column in motif_pfm_matrix.columns:
        max_disrupt_score=0
        max_R=''
        max_A=''
        indexlen=len(motif_pfm_matrix.index)
        for i in range(0,indexlen):
            for j in range(i,indexlen):
                R=motif_pfm_matrix.index[i]
                A=motif_pfm_matrix.index[j]
                scoreR=motif_pfm_matrix.loc[R,column]
                scoreA=motif_pfm_matrix.loc[A,column]
                cur_disrupt_score=scoreA-scoreR
                #print(str(cur_disrupt_score))
                #print(R)
                #print(A)
                tuple_key=(str(column),str(R),str(A))
                disrupt_dict[tuple_key]=cur_disrupt_score
                #if (scoreR>0 or scoreA>0) and abs(cur_disrupt_score)>abs(max_disrupt_score):
                #    max_disrupt_score=cur_disrupt_score
                #    max_R=R
                #    max_A=A
        #print(max_R+":::"+max_A)    
        #max_disrupt_score_list.append(max_disrupt_score)
    return disrupt_dict

def calculate_theoretic_max(motif_pfm_filename, pseudo_count=1E-4, background=None):
    motif_pfm_matrix=pd.read_csv(motif_pfm_filename, header=None,delim_whitespace=True)
    motif_pfm_matrix.index=['A','C','G','T']
    motif_pfm_matrix=motif_pfm_matrix.astype('float64')
    bg = {'A':0.2952,'C':0.2048,'G':0.2048,'T':0.2952}
    for key in bg.keys():
        bg[key]=np.log(bg[key])
    motif_pfm_matrix=motif_pfm_matrix+pseudo_count
    pfm_min=np.min(motif_pfm_matrix,axis=0)
    pfm_max=np.max(motif_pfm_matrix,axis=0)
    pfm_sum=np.sum(motif_pfm_matrix,axis=0)
    min_seq=np.log(pfm_min/pfm_sum)
    max_seq=np.log(pfm_max/pfm_sum)
    pfm_max_id=motif_pfm_matrix.idxmax()
    pfm_min_id=motif_pfm_matrix.idxmin()
    bg_min_seq=[bg[nucid] for nucid in pfm_min_id]
    bg_max_seq=[bg[nucid] for nucid in pfm_max_id]
    min_score=sum(min_seq)-sum(bg_min_seq)
    max_score=sum(max_seq)-sum(bg_max_seq)
    motif_pfm_matrix_2=np.log(motif_pfm_matrix/pfm_sum)
    max_disrupt_score_list=[]
    for column in motif_pfm_matrix_2.columns:
        max_disrupt_score=0
        max_R=''
        max_A=''
        indexlen=len(motif_pfm_matrix_2.index)
        for i in range(0,indexlen):
            for j in range(i,indexlen):
                R=motif_pfm_matrix_2.index[i]
                A=motif_pfm_matrix_2.index[j]
                scoreR=motif_pfm_matrix_2.loc[R,column]-bg[R]
                scoreA=motif_pfm_matrix_2.loc[A,column]-bg[A]
                cur_disrupt_score=scoreR-scoreA
                #print(str(cur_disrupt_score))
                #print(R)
                #print(A)
                if (scoreR>0 or scoreA>0) and abs(cur_disrupt_score)>abs(max_disrupt_score):
                    max_disrupt_score=cur_disrupt_score
                    max_R=R
                    max_A=A
        #print(max_R+":::"+max_A)    
        max_disrupt_score_list.append(max_disrupt_score)
    max_disrupt=max(np.abs(max_disrupt_score_list))
    return max_score, max_disrupt

def write_scale_file(motif_pfm_folder, motif_scan_folder, motif_id):
    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_score_filename=os.path.join(motif_scan_folder, motif_id+".scores")
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    scale_file=open(motif_scale_filename,"w")
    score_df=pd.read_csv(motif_score_filename,sep="\t")
    score_df_sub=score_df[(score_df['binding_ref']>0)&(score_df['binding_alt']>0)]
    max_binding_score=max(max(score_df['binding_ref']),max(score_df['binding_alt']))
    max_disrupt_score=max(score_df['disrupt_score'])
    theoretic_max_bind, theoretic_max_disrupt = calculate_theoretic_max(motif_pfm_filename)
    headstring="bg_max_bind_score"+"\t"+"bg_max_disrupt_score"+"\t"+"theoretic_max_bind_score"+"\t"+"theoretic_max_disrupt_score"
    scale_file.write(headstring+"\n")
    outstring=str(max_binding_score)+"\t"+str(max_disrupt_score)+"\t"+str(theoretic_max_bind)+"\t"+str(theoretic_max_disrupt)
    scale_file.write(outstring+"\n")
    scale_file.close()


def run_scale(motif_pfm_folder, motif_scan_folder, numberofthreads):
    ext = [".scores"]
    motif_id_list = []
    for file in os.listdir("./motifscanfiles"):
        if file.endswith(tuple(ext)):
            motif_id=os.path.basename(file).replace(".scores","")
            motif_id_list.append(motif_id)
    p=multiprocessing.Pool(processes=numberofthreads)
    combine_func=partial(write_scale_file,motif_pfm_folder, motif_scan_folder)
    p.map(combine_func, motif_id_list)
    p.close()


def motif_scan(sequence, motif_pfm_matrix, motif_len, pos_start, pos_stop, strand, output_dict):
    if output_dict==None:
        output_dict=dict()
    #motif_len=motif_pfm_matrix.shape[1]
    seq_len=len(sequence)
    sequence=sequence.upper()
    if strand=="+":
        for start_pos in range(pos_start, pos_stop+1):
            if start_pos+motif_len>seq_len:
                break
            this_sequence=sequence[start_pos:start_pos+motif_len]
            score=sequence_motif_score(this_sequence, motif_pfm_matrix)
            output_dict[start_pos+1]=(start_pos,strand,score)
    elif strand=="-":
        old_chars = "ACGTYSMBDRWKVHN"
        replace_chars = "TGCARWKVHYSMBDN"
        reversecomp_tab = str.maketrans(old_chars,replace_chars)
        sequence_r = sequence.translate(reversecomp_tab)[::-1]
        trans_pos_stop=seq_len-1-(pos_start+motif_len-1)
        trans_pos_start=seq_len-1-(pos_stop+motif_len-1)
        
        for start_pos in range(trans_pos_start, trans_pos_stop+1):
            if start_pos+motif_len>seq_len:
                break
            this_sequence=sequence_r[start_pos:start_pos+motif_len]
            score=sequence_motif_score(this_sequence, motif_pfm_matrix)
            trans_start_pos=seq_len-1-(start_pos+motif_len-1)
            output_dict[-(trans_start_pos+1)]=(trans_start_pos,strand,score)
        
    return output_dict


def extract_all_motifs_triplets_new(motif_pfm_matrix, motif_len, sequence, snp_start_pos,snp_end_pos):
    min_start=max(min(snp_start_pos-motif_len+1,snp_end_pos-motif_len+1),0)
    max_start=min(max(snp_start_pos,snp_end_pos),len(sequence)-1)
    output_dict=dict()
    output_dict=motif_scan(sequence, motif_pfm_matrix, motif_len, min_start, max_start, "+", output_dict)
    output_dict=motif_scan(sequence, motif_pfm_matrix, motif_len, min_start, max_start, "-", output_dict)
    
    return output_dict
def calculate_raw_area_score(binding_score, disrupt_score):
    area_score=0
    binding_score=float(binding_score)
    disrupt_score=float(disrupt_score)
    sign=np.sign(disrupt_score)
    area_score=sign*abs(binding_score*disrupt_score)
    return area_score

def calculate_area_on_the_fly_new(currentSNP, motif_pfm_folder, motif_id):
    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_pfm_matrix=generate_pfm_matrix(motif_pfm_filename, 1E-4, None)
    motif_len=motif_pfm_matrix.shape[1]
    motif_hits_ref = extract_all_motifs_triplets_new(motif_pfm_matrix, motif_len, currentSNP.seq_ref, currentSNP.pos_start_ref,currentSNP.pos_end_ref)
    motif_hits_alt = extract_all_motifs_triplets_new(motif_pfm_matrix, motif_len, currentSNP.seq_alt, currentSNP.pos_start_alt,currentSNP.pos_end_alt)
    disruption_area_score=0
    disruption_score=0
    max_binding_score=-999999
    disruption_part=(0,0,0)
    alternative_part=(0,0,0)
    offset=currentSNP.pos_end_alt-currentSNP.pos_start_alt-(currentSNP.pos_end_ref-currentSNP.pos_start_ref)
    for key_pos in motif_hits_ref.keys():
        this_part=motif_hits_ref[key_pos]
        counter_parts=[]
        conservative_one=(0,0,0)
        this_disrupt_score=0
        if key_pos>0:
            if key_pos<currentSNP.pos_start_ref:
                counter_key_pos=key_pos
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
            if key_pos+motif_len-1>currentSNP.pos_end_ref:
                counter_key_pos=key_pos+offset
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
        elif key_pos<0:
            if key_pos>currentSNP.pos_start_ref-len(currentSNP.seq_ref)+1:
                counter_key_pos=key_pos
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
            if key_pos-motif_len+1<currentSNP.pos_end_ref-len(currentSNP.seq_ref)+1:
                counter_key_pos=key_pos-offset
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])

        minimal_area_score=999999
        for counter_part in counter_parts:
            this_part_2=this_part[2]
            counter_part_2=counter_part[2]
            #if this_part[2]<0:
            #    this_part_2=0
            #if counter_part[2]<0:
            #    counter_part_2=0
            cur_max_bind=max(this_part_2,counter_part_2)
            #if this_part_2>=0 or counter_part_2>=0:
            if cur_max_bind>max_binding_score and cur_max_bind>0: #promote maximum binding for each sequence
                max_binding_score=cur_max_bind
                diff_score=counter_part_2-this_part_2
                this_area_score=calculate_raw_area_score(max(counter_part_2,this_part_2), diff_score)
                #if abs(this_area_score)<abs(minimal_area_score):
                if True:
                    minimal_area_score=this_area_score
                    conservative_one=counter_part
                    this_disrupt_score=diff_score
        if minimal_area_score!=999999: #and abs(minimal_area_score)>abs(disruption_area_score):
            disruption_area_score=minimal_area_score
            disruption_score=this_disrupt_score
            disruption_part=this_part
            alternative_part=conservative_one
            

    return disruption_area_score, disruption_score, disruption_part, alternative_part

def getSNPMotif_new(df_snps_valid, background_snp_id_list, df_motif_express, motif_pfm_folder, motif_scan_folder, motif_id):
    
    usefulcolumns=['rsid:motif', 'ID','UID','motif','expression','binding_score','disrupt_score','area_score']
    result_df = pd.DataFrame([], index=df_snps_valid.index, columns=usefulcolumns)
    
    result_df['motif']=[motif_id]*len(df_snps_valid.index)
    
    expression_level=0
    if motif_id in df_motif_express.index:
        expression_level=df_motif_express.loc[motif_id]['FPKM']
    result_df['expression']=[expression_level]*len(df_snps_valid.index)
    
    rsid_motif_list=[]
    ID_list=[]
    UID_list=[]
    binding_score_list=[]
    disrupt_score_list=[]
    area_score_list=[]
    
    for valid_index in df_snps_valid.index:
        rsid=df_snps_valid.loc[valid_index]['ID']
        uid=df_snps_valid.loc[valid_index]['UID']
        ID_list.append(rsid)
        UID_list.append(uid)
        
        this_disrupt_score=0
        rsid_motif=rsid+":"+motif_id
        rsid_motif_list.append(rsid_motif)

        currentsnp=df_snps_valid.loc[valid_index]      

        disruption_area_score, disruption_score, ref_part, alt_part= calculate_area_on_the_fly_new(currentsnp, motif_pfm_folder, motif_id)
        if ref_part==(0,0,0) or alt_part==(0,0,0):
            ref_part=(0,0,-1.1111)
            alt_part=(0,0,-1.1111)
        #if float(ref_part[2])<0 and float(alt_part[2])<0:
        #    continue
            
        #motif_best_score, motif_worst_score, shuffle_worst_score, permute_f=calculate_permutation_frequency(currentsnp,motif_pfm_filename)


        ref_part_2=float(ref_part[2])
        alt_part_2=float(alt_part[2])
        good_bind_score=max(ref_part_2,alt_part_2)

        binding_score_list.append(good_bind_score)
        disrupt_score_list.append(disruption_score)
        area_score_list.append(disruption_area_score)
        
    result_df['rsid:motif'] = rsid_motif_list
    result_df['ID'] = ID_list
    result_df['UID'] = UID_list
    result_df['binding_score'] = binding_score_list
    result_df['disrupt_score'] = disrupt_score_list
    result_df['area_score'] = area_score_list
    result_df_sub=result_df.loc[result_df['binding_score']>=0].copy()
    
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    max_disrupt=max(parameters.iloc[0]['bg_max_disrupt_score'],parameters.iloc[0]['theoretic_max_disrupt_score'])
    
    scaled_binding_score_list=np.minimum(result_df_sub['binding_score']/max_bind,[1]*len(result_df_sub.index))
    result_df_sub['scaled_binding_score']=scaled_binding_score_list
    scaled_disrupt_score_list=np.sign(result_df_sub['disrupt_score'])*np.minimum(np.abs(result_df_sub['disrupt_score']/max_disrupt),[1]*len(result_df_sub.index))
    result_df_sub['scaled_disrupt_score']=scaled_disrupt_score_list
    result_df_sub['scaled_area_score']=result_df_sub['scaled_binding_score']*result_df_sub['scaled_disrupt_score']
    result_df_sub.set_index('rsid:motif', inplace=True)
    return result_df_sub

def run_snp_motif_matching(target_snp_filename, background_snp_filename, motif_list_filename, motif_expression_filename, motif_pfm_folder, motif_scan_folder, output_filename, num_of_threads):
    #test
    #target_snp_filename="guillaumetop_snps_seq_pickle.df.txt"
    df_snps_valid=pd.read_csv(target_snp_filename, sep='\t')
    df_snps_valid['UID']=df_snps_valid['ID']+":"+df_snps_valid['REF']+":"+df_snps_valid['ALT']
    #df_snps_valid.set_index('ID', inplace=True)
    
    background_snp_id_list=None
    #if background_snp_filename!="genome":
    #    df_background=pd.read_csv(background_snp_filename, sep='\t')
    #    background_snp_id_list=df_background['ID'].tolist()
    
    motif_id_list=list()
    if motif_list_filename=="all":
        ext = [".scores"]
        motif_id_list = []
        for file in os.listdir(motif_scan_folder):
            if file.endswith(tuple(ext)):
                motif_id=os.path.basename(file).replace(".scores","")
                motif_id_list.append(motif_id)
    else:
        motif_id_pd=pd.read_csv(motif_list_filename, sep='\t', header=None)
        motif_id_list=motif_id_pd[0].tolist()

    
    #test
    #numberofthreads=2
    
    #test
    #outputfilename="result_df.txt"
    
    df_motif_express=pd.read_csv(motif_expression_filename, sep='\t')
    df_motif_express.set_index('motif', inplace=True)
    
    #result_df=getSNPMotif_new(df_snps_valid, background_snp_id_list, df_motif_express, motif_pfm_folder, motif_scan_folder, motif_id_list[0])


    
    p=multiprocessing.Pool(processes=num_of_threads)
    getSNPMotif_func=partial(getSNPMotif_new, df_snps_valid, background_snp_id_list, df_motif_express, motif_pfm_folder, motif_scan_folder)
    pool_results=p.map(getSNPMotif_func, motif_id_list)
    p.close()
    p.join()
    # merging parts processed by different processes
    result_df = pd.concat(pool_results, axis=0)
    
    result_df.to_csv(output_filename,sep='\t',index=True)


def calculate_bind_on_the_fly_new(currentSNP, motif_pfm_folder, motif_id):
    return_tuple_list=[]
    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_pfm_matrix=generate_pfm_matrix(motif_pfm_filename, 1E-4, None)
    motif_len=motif_pfm_matrix.shape[1]
    motif_hits_ref = extract_all_motifs_triplets_new(motif_pfm_matrix, motif_len, currentSNP.seq_ref, currentSNP.pos_start_ref,currentSNP.pos_end_ref)
    motif_hits_alt = extract_all_motifs_triplets_new(motif_pfm_matrix, motif_len, currentSNP.seq_alt, currentSNP.pos_start_alt,currentSNP.pos_end_alt)
    disruption_area_score=0
    disruption_score=0
    max_binding_score=-999999
    disruption_part=(0,0,0)
    alternative_part=(0,0,0)
    offset=currentSNP.pos_end_alt-currentSNP.pos_start_alt-(currentSNP.pos_end_ref-currentSNP.pos_start_ref)
    for key_pos in motif_hits_ref.keys():
        this_part=motif_hits_ref[key_pos]
        counter_parts=[]
        conservative_one=(0,0,0)
        this_disrupt_score=0
        if key_pos>0:
            if key_pos<currentSNP.pos_start_ref:
                counter_key_pos=key_pos
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
            if key_pos+motif_len-1>currentSNP.pos_end_ref:
                counter_key_pos=key_pos+offset
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
        elif key_pos<0:
            if key_pos>currentSNP.pos_start_ref-len(currentSNP.seq_ref)+1:
                counter_key_pos=key_pos
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])
            if key_pos-motif_len+1<currentSNP.pos_end_ref-len(currentSNP.seq_ref)+1:
                counter_key_pos=key_pos-offset
                if counter_key_pos in motif_hits_alt:
                    counter_parts.append(motif_hits_alt[counter_key_pos])

        minimal_area_score=999999
        for counter_part in counter_parts:
            this_part_2=this_part[2]
            counter_part_2=counter_part[2]
            #if this_part[2]<0:
            #    this_part_2=0
            #if counter_part[2]<0:
            #    counter_part_2=0
            cur_max_bind=max(this_part_2,counter_part_2)
            #if this_part_2>=0 or counter_part_2>=0:
            if cur_max_bind>max_binding_score and cur_max_bind>0: #promote maximum binding for each sequence
                max_binding_score=cur_max_bind
                diff_score=counter_part_2-this_part_2
                this_area_score=calculate_raw_area_score(max(counter_part_2,this_part_2), diff_score)
                #if abs(this_area_score)<abs(minimal_area_score):
                if True:
                    minimal_area_score=this_area_score
                    conservative_one=counter_part
                    this_disrupt_score=diff_score
        if minimal_area_score!=999999: #and abs(minimal_area_score)>abs(disruption_area_score):
            disruption_area_score=minimal_area_score
            disruption_score=this_disrupt_score
            disruption_part=this_part
            alternative_part=conservative_one
            new_bind_tuple=(this_part[0],this_part[1],max(this_part[2],conservative_one[2]))
            return_tuple_list.append(new_bind_tuple)

    return return_tuple_list

def scan_overlap_all_2(df_snps_valid, df_motif_express, motif_pfm_folder, motif_scan_folder, motif_id):
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    max_disrupt=max(parameters.iloc[0]['bg_max_disrupt_score'],parameters.iloc[0]['theoretic_max_disrupt_score'])

    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_pfm_matrix=generate_pfm_matrix(motif_pfm_filename, 1E-4, None)
    motif_len=motif_pfm_matrix.shape[1]
    overlap=dict()
    overlap['#CHR']=list()
    overlap['Start']=list()
    overlap['End']=list()
    overlap['Motif']=list()
    overlap['UID']=list()
    overlap['Score']=list()
    overlap['Strand']=list()
    overlap['Ref_score']=list()
    overlap['Alt_score']=list()
    overlap['Disruption_score']=list()
    overlap['Expression']=list()

    expression_level=0
    if motif_id in df_motif_express.index:
        expression_level=df_motif_express.loc[motif_id]['FPKM']

    for valid_index in df_snps_valid.index:
        currentSNP=df_snps_valid.loc[valid_index]
        a=int(currentSNP['pos_start_ref'])
        b=int(currentSNP['pos_end_ref'])
        left_end=min(a,b)
        #overlap_tuple_list=calculate_bind_on_the_fly_new(currentSNP, motif_pfm_folder, motif_id)

        disruption_area_score, disruption_score, ref_part, alt_part= calculate_area_on_the_fly_new(currentSNP, motif_pfm_folder, motif_id)
        if ref_part==(0,0,0) or alt_part==(0,0,0):
            continue
        else:
        #for value_tuple in overlap_tuple_list:
            value_tuple=None
            if float(ref_part[2])>float(alt_part[2]):
                value_tuple=ref_part
            else:
                value_tuple=alt_part
            bind_score=float(value_tuple[2])
            if bind_score>0:
                bind_score=min(bind_score/max_bind, 1)
                overlap['Score'].append(bind_score)
                overlap['Strand'].append(value_tuple[1])

                overlap['#CHR'].append("chr"+str(currentSNP['CHR']))
                start_pos=currentSNP['POS']-1-left_end+int(value_tuple[0])
                end_pos=start_pos+motif_len
                overlap['Start'].append(start_pos)
                overlap['End'].append(end_pos)
                overlap['Motif'].append(motif_id)
                overlap['UID'].append(currentSNP['UID'])
                overlap['Ref_score'].append(np.sign(ref_part[2])*min(np.abs(float(ref_part[2]))/max_bind,1))
                overlap['Alt_score'].append(np.sign(alt_part[2])*min(np.abs(float(alt_part[2]))/max_bind,1))
                overlap['Disruption_score'].append(np.sign(disruption_score)*np.minimum(np.abs(disruption_score/max_disrupt),1))
                overlap['Expression'].append(expression_level)
    overlap_df=pd.DataFrame(overlap)
    return overlap_df




def scan_overlap_all(df_snps_valid, motif_pfm_folder, motif_scan_folder, motif_id):
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_pfm_matrix=generate_pfm_matrix(motif_pfm_filename, 1E-4, None)
    motif_len=motif_pfm_matrix.shape[1]
    overlap=dict()
    overlap['CHR']=list()
    overlap['Start']=list()
    overlap['End']=list()
    overlap['Motif']=list()
    overlap['UID']=list()
    overlap['Score']=list()
    overlap['Strand']=list()
    for valid_index in df_snps_valid.index:
        currentSNP=df_snps_valid.loc[valid_index]
        a=int(currentSNP['pos_start_ref'])
        b=int(currentSNP['pos_end_ref'])
        left_end=min(a,b)
        overlap_tuple_list=calculate_bind_on_the_fly_new(currentSNP, motif_pfm_folder, motif_id)
        for value_tuple in overlap_tuple_list:
            bind_score=float(value_tuple[2])
            if bind_score>0:
                bind_score=min(bind_score/max_bind, 1)
                overlap['Score'].append(bind_score)
                overlap['Strand'].append(value_tuple[1])
                
                overlap['CHR'].append("chr"+str(currentSNP['CHR']))
                start_pos=currentSNP['POS']-1-left_end+int(value_tuple[0])
                end_pos=start_pos+motif_len
                overlap['Start'].append(start_pos)
                overlap['End'].append(end_pos)
                overlap['Motif'].append(motif_id)
                overlap['UID'].append(currentSNP['UID'])
    overlap_df=pd.DataFrame(overlap)
    return overlap_df


def scan_nonoverlap(left_sequence, right_sequence, motif_pfm_matrix):
    motif_len=motif_pfm_matrix.shape[1]
    min_start=0
    max_start=len(left_sequence)-motif_len
    output_dict=dict()
    output_dict=motif_scan(left_sequence, motif_pfm_matrix, motif_len, min_start, max_start, "+", output_dict)
    output_dict=motif_scan(left_sequence, motif_pfm_matrix, motif_len, min_start, max_start, "-", output_dict)
    max_start=len(right_sequence)-motif_len
    output_dict2=dict()
    output_dict2=motif_scan(right_sequence, motif_pfm_matrix, motif_len, min_start, max_start, "+", output_dict2)
    output_dict2=motif_scan(right_sequence, motif_pfm_matrix, motif_len, min_start, max_start, "-", output_dict2)
    return output_dict, output_dict2
def scan_nonoverlap_all(df_snps_valid, df_motif_express, motif_pfm_folder, motif_scan_folder, motif_id):
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    motif_pfm_filename=os.path.join(motif_pfm_folder, motif_id+".pfm")
    motif_pfm_matrix=generate_pfm_matrix(motif_pfm_filename, 1E-4, None)
    motif_len=motif_pfm_matrix.shape[1]
    non_overlap=dict()
    non_overlap['#CHR']=list()
    non_overlap['Start']=list()
    non_overlap['End']=list()
    non_overlap['Motif']=list()
    non_overlap['UID']=list()
    non_overlap['Distance']=list()
    non_overlap['Direction']=list()
    non_overlap['Score']=list()
    non_overlap['Strand']=list()
    non_overlap['Expression']=list()

    expression_level=0
    if motif_id in df_motif_express.index:
        expression_level=df_motif_express.loc[motif_id]['FPKM']

    for valid_index in df_snps_valid.index:
        currentSNP=df_snps_valid.loc[valid_index]
        sequence=currentSNP['seq_ref'].upper()
        a=int(currentSNP['pos_start_ref'])
        b=int(currentSNP['pos_end_ref'])
        left_end=min(a,b)
        right_start=max(a,b)
        left_sequence=sequence[0:left_end]
        right_sequence=sequence[right_start+1:]
        left_dict, right_dict = scan_nonoverlap(left_sequence, right_sequence, motif_pfm_matrix)
        for key in left_dict.keys():
            value_tuple=left_dict[key]
            bind_score=float(value_tuple[2])
            if bind_score>0:
                bind_score=min(bind_score/max_bind, 1)
                non_overlap['Score'].append(bind_score)
                non_overlap['Strand'].append(value_tuple[1])
                
                non_overlap['#CHR'].append("chr"+str(currentSNP['CHR']))
                distance=left_end-int(value_tuple[0])
                non_overlap['Distance'].append(distance)
                non_overlap['Direction'].append('left')
                start_pos=currentSNP['POS']-distance-1
                end_pos=start_pos+motif_len
                non_overlap['Start'].append(start_pos)
                non_overlap['End'].append(end_pos)
                non_overlap['Motif'].append(motif_id)
                non_overlap['UID'].append(currentSNP['UID'])
                non_overlap['Expression'].append(expression_level)
        for key in right_dict.keys():
            value_tuple=right_dict[key]
            bind_score=float(value_tuple[2])
            if bind_score>0:
                bind_score=min(bind_score/max_bind, 1)
                non_overlap['Score'].append(bind_score)
                non_overlap['Strand'].append(value_tuple[1])
                
                non_overlap['#CHR'].append("chr"+str(currentSNP['CHR']))
                distance=int(value_tuple[0])+1
                non_overlap['Distance'].append(distance)
                non_overlap['Direction'].append('right')
                start_pos=currentSNP['POS']+distance-1
                end_pos=start_pos+motif_len
                non_overlap['Start'].append(start_pos)
                non_overlap['End'].append(end_pos)
                non_overlap['Motif'].append(motif_id)
                non_overlap['UID'].append(currentSNP['UID'])
                non_overlap['Expression'].append(expression_level)
    non_overlap_df=pd.DataFrame(non_overlap)
    return non_overlap_df

def run_snp_motif_binding(target_snp_filename,motif_list_filename,motif_expression_filename,motif_pfm_folder,motif_scan_folder,output_dir,num_of_threads):
   
    #target_snp_filename="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_04_GeneLinker/2018_08_06/hit_snps_seq_50_pickle.df.txt"
    #motif_list_filename="test_motif_list.txt"
    #motif_expression_filename="HUDEP2_median_expression.csv"
    #motif_pfm_folder="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_08_new_data_structure/JASPAR/pfmfiles"
    #motif_scan_folder="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_08_new_data_structure/JASPAR/motifscanfiles"
    #output_dir="/data/pinello/PEOPLE/qiuming/2018_01_MOTIF_RAPTOR/2018_04_GeneLinker/2018_09_22/nonoverlap_out"
    
    df_snps_valid=pd.read_csv(target_snp_filename, sep='\t')
    df_snps_valid['UID']=df_snps_valid['ID']+":"+df_snps_valid['REF']+":"+df_snps_valid['ALT']
    
    ####change SNP file to bed file###

    usefulcolumns=['CHR', 'Start','End','UID']
    result_df = pd.DataFrame([], index=df_snps_valid['ID'], columns=usefulcolumns)
    result_df['UID']=[str(i) for i in df_snps_valid['UID']]
    chr_list=["chr"+str(i) for i in df_snps_valid['CHR']]
    start_list=[str(int(i)-1) for i in df_snps_valid['POS']]
    end_list=[str(i) for i in df_snps_valid['POS']]
    result_df['CHR']=chr_list
    result_df['Start']=start_list
    result_df['End']=end_list
    result_df.to_csv(os.path.join(output_dir,"SNP.bed"),index=None, header=None,sep='\t')



    motif_id_list=list()
    if motif_list_filename=="all":
        ext = [".scores"]
        motif_id_list = []
        for file in os.listdir(motif_scan_folder):
            if file.endswith(tuple(ext)):
                motif_id=os.path.basename(file).replace(".scores","")
                motif_id_list.append(motif_id)
    else:
        motif_id_pd=pd.read_csv(motif_list_filename, sep='\t', header=None)
        motif_id_list=motif_id_pd[0].tolist()

    nonoverlap_filename=os.path.join(output_dir,"non_overlap.bed")
    overlap_filename=os.path.join(output_dir,"overlap.bed")
    
    #motif_id=motif_id_list[0]
    
    
    #non_overlap_df=scan_nonoverlap_all(df_snps_valid, motif_pfm_folder, motif_scan_folder, motif_id)
    #overlap_df= scan_overlap_all(df_snps_valid, motif_pfm_folder, motif_scan_folder, motif_id)
    
    df_motif_express=pd.read_csv(motif_expression_filename, sep='\t')
    df_motif_express.set_index('motif', inplace=True)
    
    p=multiprocessing.Pool(processes=num_of_threads)
    getnonoverlap_func=partial(scan_nonoverlap_all, df_snps_valid, df_motif_express, motif_pfm_folder, motif_scan_folder)
    nonoverlap_pool_results=p.map(getnonoverlap_func, motif_id_list)
    p.close()
    p.join()
    # merging parts processed by different processes
    non_overlap_df = pd.concat(nonoverlap_pool_results, axis=0)
    
    non_overlap_df.to_csv(nonoverlap_filename,index=None, header=True,sep='\t')

    p=multiprocessing.Pool(processes=num_of_threads)
    getoverlap_func=partial(scan_overlap_all_2, df_snps_valid, df_motif_express, motif_pfm_folder, motif_scan_folder)
    overlap_pool_results=p.map(getoverlap_func, motif_id_list)
    p.close()
    p.join()
    # merging parts processed by different processes
    overlap_df = pd.concat(overlap_pool_results, axis=0)
    
    overlap_df.to_csv(overlap_filename,index=None, header=True,sep='\t')


def main():
    parser = argparse.ArgumentParser(prog='snp_motif_scan', description='Scan motif on SNPs on the fly.')
    subparsers = parser.add_subparsers(help='help for subcommand: scale, pfmscan, bindscan', dest="command")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    # create the parser for the "scale" command
    parser_a = subparsers.add_parser('scale', help='scale help')
    parser_a.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_a.add_argument('-score', '--score_folder', type=str, help='motif score files folder',dest="score_folder")
    parser_a.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    # create the parser for the "pfmscan" command
    parser_b = subparsers.add_parser('pfmscan', help='scan motif from pfm files')
    parser_b.add_argument('-t', '--target_snp', type=str, help='target snp list file with sequences',dest="target_snps")
    parser_b.add_argument('-bg', '--bg_snp', type=str, help='background snp list file or (genome)',dest="bg_snps")
    parser_b.add_argument('-m', '--motifs', type=str, help='motif list file, no header',dest="motif_list")
    parser_b.add_argument('-e', '--expression', type=str, help='expression file',dest="expression")
    parser_b.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_b.add_argument('-score', '--score_folder', type=str, help='motif score files folder (with scaling parameters)',dest="score_folder")
    parser_b.add_argument('-mo', '--motifscan_output', type=str, action='store', default='result_new_df.txt', help='SNPs with Motif Scan Ouput',dest="snp_motif_out")
    parser_b.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    # create the parser for the "bindscan" command
    parser_c = subparsers.add_parser('bindscan', help='only scan for binding scores, SNP overlap and nonoverlap window')
    parser_c.add_argument('-t', '--target_snp', type=str, help='target snp list file with sequences',dest="target_snps")
    parser_c.add_argument('-m', '--motifs', type=str, help='motif list file, no header',dest="motif_list")
    parser_c.add_argument('-e', '--expression', type=str, help='expression file',dest="expression")
    parser_c.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_c.add_argument('-score', '--score_folder', type=str, help='motif score files folder (with scaling parameters)',dest="score_folder")
    parser_c.add_argument('-od', '--bed_output_folder', type=str, help='Ouput Folder for bed files',dest="dir_out")
    parser_c.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")


    args = parser.parse_args()
    
    #parser.print_help()
    
    if args.command=="scale":
        print("Command: generating scale files...")
        run_scale(args.pfm_folder, args.score_folder, args.thread_num)
    elif args.command=="pfmscan":
        print("Command: scanning...")
        run_snp_motif_matching(args.target_snps,args.bg_snps,args.motif_list,args.expression,args.pfm_folder,args.score_folder,args.snp_motif_out,args.thread_num)
    elif args.command=="bindscan":
        print("Command: scanning...")
        run_snp_motif_binding(args.target_snps,args.motif_list,args.expression,args.pfm_folder,args.score_folder,args.dir_out,args.thread_num)
   


       

if __name__=="__main__":
    #sys.argv = "snp_motif_scan.py scale -pfm pfmfiles -score motifscanfiles -p 2"
    main()
