import sys, traceback
import pandas as pd
import numpy as np
import time, datetime
import matplotlib.pyplot as plt
import re
import argparse
import os
import multiprocessing
from functools import partial
from motif_matching_lcp import motif_matching

def index_func(command_address, outputdirname, seqfilename):
    samplename = os.path.basename(seqfilename).replace("_SEQ_REF_ALT.txt","")
    seq_filename = seqfilename  # This should contain the sequence to scan
    sa_filename = os.path.join(outputdirname, samplename + ".sa")   # This will contain the suffix array
    lcp_filename = os.path.join(outputdirname, samplename + ".lcp")  # This will contain the lcp array
    mksary_command=command_address+" "+seq_filename+" "+sa_filename+" "+lcp_filename
    status=os.system(mksary_command) #./mksary $seq_file $sa_file $lcp_file

def run_mksary(mksaray_path, seq_file_folder, outputdirname, numberofthreads):
    valid = re.compile(".+_SEQ_REF_ALT.txt")
    files_to_enc = []
    for root, dirs, files in os.walk(seq_file_folder):
        for file in files:
            if valid.match(file) is not None:
                files_to_enc.append(os.path.join(root, file))
    '''
    for filename in files_to_enc:
        command_address = os.path.join(mksaray_path, "mksary") # This is the command path
        samplename = os.path.basename(filename).replace("_SEQ_REF_ALT.txt","")
        seq_filename = filename  # This should contain the sequence to scan
        sa_filename = os.path.join(outputdirname, samplename + ".sa")   # This will contain the suffix array
        lcp_filename = os.path.join(outputdirname, samplename + ".lcp")  # This will contain the lcp array
        mksary_command=command_address+" "+seq_filename+" "+sa_filename+" "+lcp_filename
        status=os.system(mksary_command) #./mksary $seq_file $sa_file $lcp_file
    '''
    command_address = os.path.join(mksaray_path, "mksary") # This is the command path

    p=multiprocessing.Pool(processes=numberofthreads)
    whole_index_func=partial(index_func, command_address, outputdirname)
    p.map(whole_index_func, files_to_enc)
    p.close()

def generate_pfm_matrix(pseudo_count, background, motif_pfm_filename):
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
    
    
    return motif_pfm_matrix_2.values

def load_files(seq_filename, sa_filename, lcp_filename):
    # Load the sequence text file (assumes upper case letters)
    #print("Start file loading...  ")
    start_time = time.time()
    Ftext = open(seq_filename, "r")
    text = Ftext.read()
    n = len(text)
    #print(str(n) + "bases, loaded in "+ str((time.time() - start_time)) + " seconds.")


    # Load the SA
    #print("Start SA loading...  ")
    start_time = time.time()
    Fsa = open(sa_filename, "rb")
    dt = np.dtype('u4')
    sa = np.fromfile(Fsa, dtype=dt, count = n, sep ="")
    #print("in " + str((time.time() - start_time)) + " seconds." )

    # Load the LCP
    #print("Start LCP loading...  ")
    start_time = time.time()
    Flcp = open(lcp_filename, "rb")
    dt = np.dtype('u4')
    lcp = np.fromfile(Flcp, dtype=dt, count = n, sep ="")
    #print("in " + str((time.time() - start_time)) + " seconds." )

    # Printing the statistics (change the boolean if you wish to see the statistics)

    if (False):
        print("Plotting statistics... ")
        new_lcp = lcp.copy()
        new_lcp.sort()
        maxlcp = np.max(new_lcp)
        #print(" max lcp = ", maxlcp)

        plt.yscale('log')
        plt.title('LCP distribution (bucket size = 10)')
        plt.xlabel('LCP length')
        plt.ylabel('frequency')
        plt.hist(new_lcp, color = 'blue', edgecolor = 'black', bins = int(maxlcp/10))
        plt.show()

    #print("Loading files done!")
    return n, text, sa, lcp

def scan_motif(n, text, sa, lcp, pos_filename, outputdirname, motif_pfm_filename_matrix_tuple):
    motif_pfm_filename = motif_pfm_filename_matrix_tuple[0]
    motif_pfm_matrix = motif_pfm_filename_matrix_tuple[1]
    block_size = 61 #size of the window centered on the SNPs
    num_tot_blocks = ( n // block_size )  # blocks of size 61 bases 
    num_blocks = num_tot_blocks // 2  # blocks of size 61 bases in each REF and ALT array
    
    samplename = os.path.basename(pos_filename).replace("_POS_REF_ALT.txt","")
    motifname = os.path.basename(motif_pfm_filename).replace(".pfm","")
    output_file_name  = motifname+"__"+samplename + ".scores"  # This will contain the output of the program
    output_file_name  = os.path.join(outputdirname, motifname, output_file_name)
    outputfolder = os.path.dirname(output_file_name)
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    length_motif, scores_in_ref, scores_in_alt = motif_matching(n, num_blocks, block_size, lcp.astype(np.uint), sa.astype(np.uint), text, motif_pfm_matrix)


    # Load the starting positions of the 61-size blocks
    #print("Loading block starting...  ")
    start_time = time.time()
    ref_starting = np.zeros((num_blocks), dtype = np.int)
    alt_starting = np.zeros((num_blocks), dtype = np.int)

    nb = 0
    with open(pos_filename, "r") as f:
        lines=f.readlines()
        for line in lines:
            x = int(re.sub("\D", "", line.strip('\n')))
            if(x):
                if (nb%2 == 0):
                    ref_starting[nb // 2] = int(x)
                else:
                    alt_starting[nb // 2] = int(x)
                nb = nb + 1


    #print("    num blocks loaded                   = "+  str(nb)) 
    #print("    num blocks expected (len file/2*61) = "+ str(num_tot_blocks))
    #print("done in " +str((time.time() - start_time))+" seconds. ")

    #print("Computing the final scores and writing them on file...  ")
    start_time = time.time()

    output_file = open(output_file_name, 'w+')
    output_file.write("Genome\tpos_SNP\tSNP_in_ref\tSNP_in_alt")
    output_file.write("\tbinding_ref\tbinding_alt\tpos_disrupt\tdisrupt_score\n")    
    #output_file.close()

    # We scan all SNPs and the computed scores in REF and ALT
    #output_file = open(output_file_name, 'a+')

    for snp in range(0,num_blocks):  

        # dummy initializations for the max calculation
        max_score_ref = float('-inf')  
        max_pos_ref = -1
        max_score_alt = float('-inf')
        max_pos_alt = -1

        snp_position = ref_starting[snp] + 30  # We assume = alt_starting[snp]+30
        char_SNP_ref = text[ (2 * snp) * 61 + 30 ]
        char_SNP_alt = text[ (2 * snp + 1) * 61 + 30 ]

        # includes the SNP position, which is at 30 (we count from 0)
        # start from a position that guarantees the overlapping between motif and SNP
        for i in range(30 - length_motif + 1, 31): 

            # Compute the max of the REF scores
            if (scores_in_ref[snp][i] > max_score_ref):
                max_score_ref = scores_in_ref[snp][i]
                max_pos_ref = i

            # Compute the max of the ALT scores
            if (scores_in_alt[snp][i] > max_score_alt):
                max_score_alt = scores_in_alt[snp][i]
                max_pos_alt = i


        if (max_score_ref >= max_score_alt): # the binding score is larger for REF
            disruption_score = max_score_ref - scores_in_alt[snp][max_pos_ref]
            disruption_pos = ref_starting[snp] + max_pos_ref   # absolute in the original sequence
            max1 = max_score_ref
            max2 = scores_in_alt[snp][max_pos_ref]

        else:   # the binding score is larger for ALT
            disruption_score = max_score_alt - scores_in_ref[snp][max_pos_alt]
            disruption_pos = alt_starting[snp] + max_pos_alt   # absolute in the original sequence
            max1 = scores_in_ref[snp][max_pos_alt]
            max2 = max_score_alt


        output_file.write(samplename + "\t" + str(snp_position) + "\t" + char_SNP_ref + "\t" + char_SNP_alt + "\t")
        output_file.write(str(max1) + "\t" + str(max2) + "\t" + str(disruption_pos) + "\t" + str(disruption_score) + "\n")    

    output_file.close()
    #print("Done scanning "+motifname+" in "+  str((time.time() - start_time)) + " seconds.")

def combine_clean_motifscanfiles(genome_database_folder, outputdirname, motifdirname):
    motifname=os.path.basename(motifdirname)
    combined_motif_score_filename=os.path.join(outputdirname, motifname+".scores")
    combined_file=open(combined_motif_score_filename,"w")
    paste_folder=os.path.join(outputdirname, motifname)
    ext = [".scores"]
    score_filelist=[]
    for root, dirs, files in os.walk(paste_folder):
        for file in files:
            if file.endswith(tuple(ext)):
                motif_score_filename = os.path.join(root, file)
                score_filelist.append(motif_score_filename)
    first_file_count=1       
    for motif_score_filename in sorted(score_filelist):
        pos1=str.rfind(motif_score_filename,"__")
        pos2=str.rfind(motif_score_filename,".scores")
        samplename=motif_score_filename[pos1+2:pos2]
        uid_filename=os.path.join(genome_database_folder, samplename+".uid")
        uid_file=open(uid_filename,"r")
        score_file=open(motif_score_filename,"r")
        score_line=score_file.readline().strip()
        uid_line=uid_file.readline().strip()
        line_count=1
        while(score_line!=''):
            if(line_count>1 or first_file_count==1):
                combined_line=score_line+"\t"+uid_line+"\n"
                combined_file.write(combined_line)
            score_line=score_file.readline().strip()
            uid_line=uid_file.readline().strip()
            line_count=line_count+1
        
        uid_file.close()
        score_file.close()
        first_file_count=first_file_count+1
        
    combined_file.close()

def run_scan_motif(genome_database_folder, motif_pfm_folder, outputdirname, numberofthreads):
    #get all motif pfm files
    ext = [".pfm"]
    #motif_files_to_enc = []
    #motif_matrix_list = []
    motif_files_matrix_tuple_list = []
    for root, dirs, files in os.walk(motif_pfm_folder):
        for file in files:
            if file.endswith(tuple(ext)):
                motif_pfm_filename = os.path.join(root, file)
                #motif_files_to_enc.append(motif_pfm_filename)
            
                pseudo_count=0.001
                background=None
                motif_pfm_matrix = generate_pfm_matrix(pseudo_count, background, motif_pfm_filename)
                #motif_matrix_list.append(motif_pfm_matrix)
                motif_files_matrix_tuple = (motif_pfm_filename, motif_pfm_matrix)
                motif_files_matrix_tuple_list.append(motif_files_matrix_tuple)
    
    #get all chromosome/sample indexed files
    ext = [".sa"]
    sa_files_to_enc = []
    for root, dirs, files in os.walk(genome_database_folder):
        for file in files:
            if file.endswith(tuple(ext)):
                sa_files_to_enc.append(os.path.join(root, file))

    for safile in sa_files_to_enc:
        samplename = os.path.basename(safile).replace(".sa","")
        seq_filename = os.path.join(genome_database_folder, samplename+"_SEQ_REF_ALT.txt")
        pos_filename = os.path.join(genome_database_folder, samplename+"_POS_REF_ALT.txt")
        sa_filename = os.path.join(genome_database_folder, samplename+".sa")
        lcp_filename = os.path.join(genome_database_folder, samplename+".lcp")
        n, text, sa, lcp = load_files(seq_filename, sa_filename, lcp_filename)
        
        
        p=multiprocessing.Pool(processes=numberofthreads)
        whole_scan_func=partial(scan_motif, n, text, sa, lcp, pos_filename, outputdirname)
        p.map(whole_scan_func, motif_files_matrix_tuple_list)
        p.close()
        
    ###combine files by uid, and for all chromosomes
    dirs_to_enc = []
    for root, dirs, files in os.walk(outputdirname):
        for dirname in dirs:
            dirs_to_enc.append(os.path.join(root, dirname))
    p=multiprocessing.Pool(processes=numberofthreads)
    combine_func=partial(combine_clean_motifscanfiles,genome_database_folder, outputdirname)
    p.map(combine_func, dirs_to_enc)
    p.close()


def main():
    parser = argparse.ArgumentParser(prog='Motif_Scan', description='Build genome index, run motif scan.')
    subparsers = parser.add_subparsers(help='help for subcommand: index, pfmscan', dest="command")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    # create the parser for the "index" command
    parser_a = subparsers.add_parser('index', help='index help')
    parser_a.add_argument('-g', '--genome_db', type=str, help='genome_database_folder',dest="genome_db")
    parser_a.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    # create the parser for the "pfmscan" command
    parser_b = subparsers.add_parser('pfmscan', help='scan motif from pfm files')
    parser_b.add_argument('-gi', '--indexed_genome_db', type=str, help='indexed genome_database_folder',dest="indexed_genome_db")
    parser_b.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_b.add_argument('-mo', '--motifscan_output', type=str, action='store', default='./motifscanfiles', help='Motif Scan Ouput Folder',dest="outputmotifscandir")
    parser_b.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")
    
    args = parser.parse_args()
    
    #parser.print_help()
    
    if args.command=="index":
        print("Command: indexing...")
        run_mksary("./", args.genome_db, args.genome_db, args.thread_num)
    elif args.command=="pfmscan":
        print("Command: scanning...")
        run_scan_motif(args.indexed_genome_db, args.pfm_folder, args.outputmotifscandir, args.thread_num)
        

if __name__=="__main__":
    start_time = time.time()
    print(str(datetime.datetime.now()))
    #sys.argv = "Motif_Scan.py index -g concatseq -p 2"
    #sys.argv = "Motif_Scan.py pfmscan -gi concatseq -pfm ./test_motifs -mo ./motifscanfiles2 -p 2"
    main()
    print(str(datetime.datetime.now()))
    print("Finished in "+  str((time.time() - start_time)/3600) + " hours.")
