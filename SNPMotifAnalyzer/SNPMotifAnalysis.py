import twobitreader
import numpy as np
import pandas as pd
import os
import argparse
import multiprocessing
from functools import partial
from scipy.stats import norm

def area_pvalue_per_motif_from_background(snp_motif_result_df, bg_SNP_df, df_motif_express, motif_scan_folder, outputdir, motif_id):
    #if 'gata' in motif_id.lower() or 'klf1' in motif_id.lower() or 'sox2' in motif_id.lower():
    usefulcolumns=['motif', 'target_num','background_num','expression','disrupt_pvalue','enhance_pvalue','abs_area_score_pvalue']
    test_result_df = pd.DataFrame([], index=[0], columns=usefulcolumns)

    #if not os.path.exists(score_filename):
    #    return test_result_df

    if not os.path.exists(outputdir):
        try:
            os.mkdir(outputdir)
        except OSError:
            pass

    output_background_folder=os.path.join(outputdir,"background_files")
    if not os.path.exists(output_background_folder):
        try:
            os.mkdir(output_background_folder)
        except OSError:
            pass

    motif_score_filename=os.path.join(motif_scan_folder, motif_id+".scores")
    SNP_score_dataframe=pd.read_csv(motif_score_filename,sep="\t")
    SNP_score_dataframe.set_index('UID', inplace=True)
    score_df=SNP_score_dataframe[(SNP_score_dataframe['binding_ref']>0)|(SNP_score_dataframe['binding_alt']>0)]

    if bg_SNP_df is not None:
        score_df_sub=score_df.merge(bg_SNP_df, on=['UID','ID'], how='inner')
    else:
        score_df_sub=score_df
    bg_score_filename=os.path.join(output_background_folder,motif_id+".scores")

    score_df_sub=score_df_sub.copy()
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    max_disrupt=max(parameters.iloc[0]['bg_max_disrupt_score'],parameters.iloc[0]['theoretic_max_disrupt_score'])
    x=np.maximum(score_df_sub['binding_alt'],score_df_sub['binding_ref'])
    y=score_df_sub['binding_alt']-score_df_sub['binding_ref']
    scaled_binding_score_list=np.minimum(x/max_bind,[1]*len(x))
    score_df_sub['scaled_binding_score']=scaled_binding_score_list
    scaled_disrupt_score_list=np.sign(y)*np.minimum(np.abs(y/max_disrupt),[1]*len(y))
    score_df_sub['scaled_disrupt_score']=scaled_disrupt_score_list
    score_df_sub['scaled_area_score']=score_df_sub['scaled_binding_score']*score_df_sub['scaled_disrupt_score']

    score_df_sub.to_csv(bg_score_filename,sep='\t',index=True)

    target_data=snp_motif_result_df.loc[snp_motif_result_df['motif']==motif_id]['scaled_area_score']
    bg_data=score_df_sub['scaled_area_score']
    #report_pvalue=subsampling(target_data,bg_data)
    pop_mean=np.mean(bg_data)
    pop_var=np.var(bg_data)
    test_mean=pop_mean

    test_var=pop_var/(len(target_data))

    if test_var==0:
        test_var=1E-20
    t=(np.mean(target_data)-test_mean)/((test_var)**(0.5))
    report_pvalue_left=norm.cdf(t)
    report_pvalue_right=1-report_pvalue_left
    #print(str(test_mean))
    #print(str(test_var)) 
    abs_target_data=abs(target_data)
    abs_bg_data=abs(bg_data)
    abs_pop_mean=np.mean(abs_bg_data)
    abs_pop_var=np.var(abs_bg_data)
    abs_test_mean=abs_pop_mean
    abs_test_var=abs_pop_var/(len(abs_target_data))
    if abs_test_var==0:
        abs_test_var=1E-20
    abs_t=(np.mean(abs_target_data)-abs_test_mean)/((abs_test_var)**(0.5))
    report_pvalue_abs=1-norm.cdf(abs_t)
    '''
    pos_target=target_data.loc[target_data>0]
    neg_target=target_data.loc[target_data<0]
    pos_bg=bg_data.loc[bg_data>0]
    neg_bg=bg_data.loc[bg_data<0]
    pos_pop_mean=np.mean(pos_bg)
    pos_pop_var=np.var(pos_bg)
    pos_test_mean=pos_pop_mean
    pos_num=len(pos_target)
    if pos_num==0:
        pos_num=pos_num+1
    pos_test_var=pos_pop_var/(pos_num)
    if pos_test_var==0:
        pos_test_var=1E-20
    pos_t=(np.mean(pos_target)-pos_test_mean)/((pos_test_var)**(0.5))
    report_pvalue_disrupt=norm.cdf(pos_t)

    neg_pop_mean=np.mean(neg_bg)
    neg_pop_var=np.var(neg_bg)
    neg_test_mean=neg_pop_mean
    neg_num=len(neg_target)
    if neg_num==0:
        neg_num=neg_num+1
    neg_test_var=neg_pop_var/(neg_num)
    if neg_test_var==0:
        neg_test_var=1E-20
    neg_t=(np.mean(neg_target)-neg_test_mean)/((neg_test_var)**(0.5))
    report_pvalue_enhance=1-norm.cdf(neg_t)
    '''

    expression_level=0
    if motif_id in df_motif_express.index:
        expression_level=df_motif_express.loc[motif_id]['FPKM']

    test_result_df['motif']=[motif_id]
    test_result_df['target_num']=[str(len(target_data))]
    test_result_df['background_num']=[str(len(bg_data))]
    test_result_df['expression']=[str(expression_level)]
    test_result_df['disrupt_pvalue']=[str(report_pvalue_left)]
    test_result_df['enhance_pvalue']=[str(report_pvalue_right)]
    test_result_df['abs_area_score_pvalue']=[str(report_pvalue_abs)]
    #test_result_df['disrupt_pvalue']=[str(report_pvalue_disrupt)]
    #test_result_df['enhance_pvalue']=[str(report_pvalue_enhance)]

    outputfilename=os.path.join(output_background_folder, motif_id+".pvalue")
    test_result_df.to_csv(outputfilename,sep='\t',index=None)
    test_result_df.set_index('motif', inplace=True)
    return test_result_df

def area_pvalue_per_motif_from_background_2(snp_motif_result_df, bg_SNP_df, df_motif_express, motif_scan_folder, outputdir, motif_id):
    #if 'gata' in motif_id.lower() or 'klf1' in motif_id.lower() or 'sox2' in motif_id.lower():
    usefulcolumns=['motif', 'target_num','background_num','expression','leftside_pvalue','rightside_pvalue','disrupt_pvalue','enhance_pvalue','abs_area_score_pvalue']
    test_result_df = pd.DataFrame([], index=[0], columns=usefulcolumns)

    #if not os.path.exists(score_filename):
    #    return test_result_df
   
    if not os.path.exists(outputdir):
        try:
            os.mkdir(outputdir)
        except OSError:
            pass
 
    output_background_folder=os.path.join(outputdir,"background_files")
    if not os.path.exists(output_background_folder):
        try:
            os.mkdir(output_background_folder)
        except OSError:
            pass    

        
    motif_score_filename=os.path.join(motif_scan_folder, motif_id+".scores")
    SNP_score_dataframe=pd.read_csv(motif_score_filename,sep="\t")
    SNP_score_dataframe.set_index('UID', inplace=True)
    score_df=SNP_score_dataframe[(SNP_score_dataframe['binding_ref']>0)|(SNP_score_dataframe['binding_alt']>0)]

    if bg_SNP_df is not None:
        score_df_sub=score_df.merge(bg_SNP_df, on=['UID','ID'], how='inner')
    else:
        score_df_sub=score_df
    bg_score_filename=os.path.join(output_background_folder,motif_id+".scores")
    
    score_df_sub=score_df_sub.copy()
    motif_scale_filename=os.path.join(motif_scan_folder, motif_id+".scale")
    parameters=pd.read_csv(motif_scale_filename, sep='\t')
    max_bind=max(parameters.iloc[0]['bg_max_bind_score'],parameters.iloc[0]['theoretic_max_bind_score'])
    max_disrupt=max(parameters.iloc[0]['bg_max_disrupt_score'],parameters.iloc[0]['theoretic_max_disrupt_score'])
    x=np.maximum(score_df_sub['binding_alt'],score_df_sub['binding_ref'])
    y=score_df_sub['binding_alt']-score_df_sub['binding_ref']
    scaled_binding_score_list=np.minimum(x/max_bind,[1]*len(x))
    score_df_sub['scaled_binding_score']=scaled_binding_score_list
    scaled_disrupt_score_list=np.sign(y)*np.minimum(np.abs(y/max_disrupt),[1]*len(y))
    score_df_sub['scaled_disrupt_score']=scaled_disrupt_score_list
    score_df_sub['scaled_area_score']=score_df_sub['scaled_binding_score']*score_df_sub['scaled_disrupt_score']

    score_df_sub.to_csv(bg_score_filename,sep='\t',index=True)

    target_data=snp_motif_result_df.loc[snp_motif_result_df['motif']==motif_id]['scaled_area_score']
    bg_data=score_df_sub['scaled_area_score']
    #report_pvalue=subsampling(target_data,bg_data)
    pop_mean=np.mean(bg_data)
    pop_var=np.var(bg_data)
    test_mean=pop_mean

    test_var=pop_var/(len(target_data))
    if test_var==0:
        test_var=1E-20
    t=(np.mean(target_data)-test_mean)/((test_var)**(0.5))
    report_pvalue_left=norm.cdf(t)
    report_pvalue_right=1-report_pvalue_left
    #print(str(test_mean))
    #print(str(test_var)) 
    abs_target_data=abs(target_data)
    abs_bg_data=abs(bg_data)
    abs_pop_mean=np.mean(abs_bg_data)
    abs_pop_var=np.var(abs_bg_data)
    abs_test_mean=abs_pop_mean
    abs_test_var=abs_pop_var/(len(abs_target_data))
    if abs_test_var==0:
        abs_test_var=1E-20
    abs_t=(np.mean(abs_target_data)-abs_test_mean)/((abs_test_var)**(0.5))
    report_pvalue_abs=1-norm.cdf(abs_t)
   
    pos_target=target_data.loc[target_data>0]
    neg_target=target_data.loc[target_data<0]
    pos_bg=bg_data.loc[bg_data>0]
    neg_bg=bg_data.loc[bg_data<0]
    pos_pop_mean=np.mean(pos_bg)
    pos_pop_var=np.var(pos_bg)
    pos_test_mean=pos_pop_mean
    pos_num=len(pos_target)
    if pos_num==0:
        pos_num=pos_num+1
    pos_test_var=pos_pop_var/(pos_num)
    if pos_test_var==0:
        pos_test_var=1E-20
    pos_t=(np.mean(pos_target)-pos_test_mean)/((pos_test_var)**(0.5))
    report_pvalue_disrupt=norm.cdf(pos_t)

    neg_pop_mean=np.mean(neg_bg)
    neg_pop_var=np.var(neg_bg)
    neg_test_mean=neg_pop_mean
    neg_num=len(neg_target)
    if neg_num==0:
        neg_num=neg_num+1
    neg_test_var=neg_pop_var/(neg_num)
    if neg_test_var==0:
        neg_test_var=1E-20
    neg_t=(np.mean(neg_target)-neg_test_mean)/((neg_test_var)**(0.5))
    report_pvalue_enhance=1-norm.cdf(neg_t)

 
    expression_level=0
    if motif_id in df_motif_express.index:
        expression_level=df_motif_express.loc[motif_id]['FPKM']
        
    test_result_df['motif']=[motif_id]
    test_result_df['target_num']=[str(len(target_data))]
    test_result_df['background_num']=[str(len(bg_data))]
    test_result_df['expression']=[str(expression_level)]
    test_result_df['leftside_pvalue']=[str(report_pvalue_left)]
    test_result_df['rightside_pvalue']=[str(report_pvalue_right)]
    test_result_df['abs_area_score_pvalue']=[str(report_pvalue_abs)]
    test_result_df['disrupt_pvalue']=[str(report_pvalue_disrupt)]
    test_result_df['enhance_pvalue']=[str(report_pvalue_enhance)]



    outputfilename=os.path.join(output_background_folder, motif_id+".pvalue")
    test_result_df.to_csv(outputfilename,sep='\t',index=None)
    test_result_df.set_index('motif', inplace=True)
    return test_result_df
                   
def run_motif_test(snp_motif_result_filename,snp_background_filename, motif_expression_filename, motif_scan_folder, outputdir, num_of_threads):
    snp_motif_result_df=pd.read_csv(snp_motif_result_filename, sep='\t')
    snp_motif_result_df.set_index('rsid:motif', inplace=True)
    bg_SNP_df=None
    if snp_background_filename!="genome":
        bg_SNP_df=pd.read_csv(snp_background_filename, sep='\t')
        bg_SNP_df['UID']=bg_SNP_df['ID']+":"+bg_SNP_df['REF']+":"+bg_SNP_df['ALT']
        bg_SNP_df.set_index('UID', inplace=True)    

    df_motif_express=pd.read_csv(motif_expression_filename, sep='\t')
    df_motif_express.set_index('motif', inplace=True)
                   
    motif_id_list=list(set(snp_motif_result_df['motif']))
    #area_pvalue_per_motif_from_background(snp_motif_result_df, bg_SNP_df, df_motif_express, motif_scan_folder, outputdir, motif_id_list[0])
    #motif_id_list=motif_id_list[0:2]
    
    p=multiprocessing.Pool(processes=num_of_threads)
    whole_create_func=partial(area_pvalue_per_motif_from_background, snp_motif_result_df, bg_SNP_df, df_motif_express, motif_scan_folder, outputdir)
    pool_results=p.map(whole_create_func, motif_id_list)
    p.close()
    p.join()
    # merging parts processed by different processes
    result_df = pd.concat(pool_results, axis=0)
    output_filename=os.path.join(outputdir,"all_motifs.pvalue")
    result_df.to_csv(output_filename,sep='\t',index=True)


def main():
    parser = argparse.ArgumentParser(prog='snp_motif_analysis', description='Analize motifs and SNPs in the dataset.')
    subparsers = parser.add_subparsers(help='help for subcommand: motif_test, snp_test', dest="command")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    # create the parser for the "scale" command
    parser_a = subparsers.add_parser('motif_test', help='motif test help')
    parser_a.add_argument('-r', '--snp_motif_result', type=str, help='snp motif result file from snp_motif_scan',dest="snp_motif_result")
    parser_a.add_argument('-bg', '--bg_snp', type=str, help='background snp list file or (genome)',dest="bg_snps")
    parser_a.add_argument('-e', '--expression', type=str, help='expression file',dest="expression")
    parser_a.add_argument('-score', '--score_folder', type=str, help='motif score files folder',dest="score_folder")
    parser_a.add_argument('-od', '--outputdir', type=str, help='Ouput folder',dest="outputdir")
    parser_a.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    # create the parser for the "pfmscan" command
    parser_b = subparsers.add_parser('snp_test', help='snp test help')
    parser_b.add_argument('-t', '--target_snp', type=str, help='target snp list file with sequences',dest="target_snps")
    parser_b.add_argument('-bg', '--bg_snp', type=str, help='background snp list file or (genome)',dest="bg_snps")
    parser_b.add_argument('-m', '--motifs', type=str, help='motif list file, no header',dest="motif_list")
    parser_b.add_argument('-e', '--expression', type=str, help='expression file',dest="expression")
    parser_b.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_b.add_argument('-score', '--score_folder', type=str, help='motif score files folder (with scaling parameters)',dest="score_folder")
    parser_b.add_argument('-mo', '--motifscan_output', type=str, action='store', default='result_new_df.txt', help='SNPs with Motif Scan Ouput',dest="snp_motif_out")
    parser_b.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")
    
    args = parser.parse_args()
    
    #parser.print_help()
    
    if args.command=="motif_test":
        print("Command: testing per motif...")
        run_motif_test(args.snp_motif_result,args.bg_snps, args.expression, args.score_folder, args.outputdir, args.thread_num)
    elif args.command=="snp_test":
        print("Command: testing per snp...")
        run_snp_motif_matching(args.target_snps,args.bg_snps,args.motif_list,args.expression,args.pfm_folder,args.score_folder,args.snp_motif_out,args.thread_num)


if __name__=="__main__":
    main()
