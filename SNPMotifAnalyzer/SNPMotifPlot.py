import numpy as np
import pandas as pd
import os
import scipy.stats as stats
import pybedtools
import re
import multiprocessing
from functools import partial
import matplotlib
import matplotlib as mpl
mpl.rc('pdf', fonttype=42)
import matplotlib.pyplot
#matplotlib.use('agg')
matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('white')
sns.set_context('poster')

from statsmodels.sandbox.stats.multicomp import multipletests
import scipy

#all_expression=pd.read_csv("all_expression.csv",sep="\t")
def calculate_expression_pvalue(all_expression, expression, motif_id):
    
    expression_list=all_expression.loc[all_expression['motif']==motif_id]
    expression_value_list0=expression_list[expression_list.columns[1:expression_list.shape[1]]].values.tolist()[0]
    expression_value_list = [x for x in expression_value_list0 if str(x) != 'nan']
    m=np.mean(expression_value_list)
    s=np.std(expression_value_list)
    if s==0:
        s=1
    z_score=(expression-m)/s
    p_value=1-scipy.stats.norm.cdf(z_score)
    return z_score, p_value

def plot_radar_df(radar_df,savefilename):
    categories=list(radar_df)[1:]
    N = len(categories)
    plt.style.use('default')
    
    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values=radar_df.iloc[0].drop('asgroup').values.flatten().tolist()
    values += values[:1]
    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)

    # Draw one axe per variable + add labels labels yet
    ticks=angles[:-1]
    plt.xticks(ticks, categories, color='black', size=10)
    for label,rot in zip(ax.get_xticklabels(),ticks):
        label.set_rotation(rot*180./np.pi)
        label.set_horizontalalignment("left")
        label.set_rotation_mode("anchor")
    
    # Draw ylabels
    ax.set_rlabel_position(15)
    #plt.yticks([10,20,30], ["10","20","30"], color="grey", size=7)
    #plt.ylim(0,40)

    # Plot data
    ax.plot(angles, values, linewidth=2, linestyle='solid')

    # Fill area
    ax.fill(angles, values, 'b', alpha=0.3)
    
    ax.grid(color='grey', linestyle='--', linewidth=1)
    plt.savefig(savefilename, bbox_inches='tight', pad_inches=0)
    #plt.show()

def plot_motif_snp_pair_main(snp_motif_result_file, all_expression_file,snp_id, motif_id_name, rsid_motifid, pdffilename, snpfeaturefile):
    snp_motif_result_df=pd.read_csv(snp_motif_result_file, sep='\t')
    hit_df_sub=snp_motif_result_df.loc[snp_motif_result_df['ID']==snp_id]

    all_expression=pd.read_csv(all_expression_file,sep="\t")
    rsid_df=hit_df_sub.loc[hit_df_sub['rsid:motif']==rsid_motifid]


    a=rsid_df['scaled_binding_score'].values[0]

    b=rsid_df['scaled_disrupt_score'].values[0]

    z_score, p_value= calculate_expression_pvalue(all_expression, rsid_df.iloc[0]['expression'], rsid_df.iloc[0]['motif'])
    d= 1-p_value

    snpfeature_df=pd.read_csv(snpfeaturefile, sep="\t", header=0)

    radar_dict={
    'asgroup': [rsid_motifid],
    'Binding': [a],
    'Disruption': [abs(b)],
    'Expression': [d]
    }

    for column_num in range(4,snpfeature_df.shape[1]):
        column_name=snpfeature_df.columns[column_num]
        score=snpfeature_df[snpfeature_df['snpid']==snp_id][column_name].values[0]
        max_score=snpfeature_df[snpfeature_df['snpid']=='max'][column_name].values[0]
        score=min(max(0,score)/max_score,1)
        radar_dict[column_name]=score



    radar_df = pd.DataFrame(radar_dict)

    plot_radar_df(radar_df,pdffilename)

'''
def plot_motif_snp_pair_main(snp_motif_result_file, all_expression_file,snp_id, motif_id_name, rsid_motifid, pdffilename, conservation_folder):
    snp_motif_result_df=pd.read_csv(snp_motif_result_file, sep='\t')
    hit_df_sub=snp_motif_result_df.loc[snp_motif_result_df['ID']==snp_id]

    all_expression=pd.read_csv(all_expression_file,sep="\t")
    rsid_df=hit_df_sub.loc[hit_df_sub['rsid:motif']==rsid_motifid]


    a=rsid_df['scaled_binding_score'].values[0]

    b=rsid_df['scaled_disrupt_score'].values[0]

    z_score, p_value= calculate_expression_pvalue(all_expression, rsid_df.iloc[0]['expression'], rsid_df.iloc[0]['motif'])
    d= 1-p_value

    uid=rsid_df['UID'].values[0]
    #c=8.97467/37.0143
    #e=0.0573333
    CADD_df=pd.read_csv(os.path.join(conservation_folder,"hit_SNP_CADD.tab"),sep="\t",header=None)
    c=0
    try:
        c=CADD_df.loc[CADD_df[0]==uid][5].values[0]
        c=c/37.0143
    except:
        pass
    phastcons_46_df=pd.read_csv(os.path.join(conservation_folder,"hit_SNP_46way.tab"),sep="\t",header=None)
    e=0
    try:
        e=phastcons_46_df.loc[phastcons_46_df[0]==uid][5].values[0]
    except:
        pass

    radar_df = pd.DataFrame({
    'asgroup': [rsid_motifid],
    'Binding': [a],
    'Disruption': [abs(b)],
    'Expression': [d],
    'Phastcons': [e],
    'CADD': [c]
    })

    plot_radar_df(radar_df,pdffilename)
'''

def plot_snp_specific_main(snp_motif_result_file, snp_id, all_expression_file,pdffilename):
    all_expression=pd.read_csv(all_expression_file,sep="\t")
    snp_motif_result_df=pd.read_csv(snp_motif_result_file, sep='\t')
    hit_df_sub=snp_motif_result_df.loc[snp_motif_result_df['ID']==snp_id]
    bg_bind_scaled_score_list=hit_df_sub['scaled_binding_score']
    bg_disrupt_scaled_score_list=hit_df_sub['scaled_disrupt_score']

    color_z_score_list=[]
    for ind in hit_df_sub.index:
        expression=hit_df_sub.loc[ind,'expression']
        motif_id=hit_df_sub.loc[ind,'motif']
        z_score, p_value=calculate_expression_pvalue(all_expression, expression, motif_id)
        color_z_score_list.append(z_score)
    #color_score=hit_df_sub['expression']
    #plt.scatter(list_b, list_a, s=4,color='blue')
    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w', edgecolor='w')
    ax = fig.add_subplot(111)
    plt.grid()
    ax.set_aspect('equal')
    #plt.scatter(bg_disrupt_scaled_score_list,bg_bind_scaled_score_list,  s=50,color='red',marker='s')
    plt.scatter(bg_disrupt_scaled_score_list,bg_bind_scaled_score_list,  s=50,c=color_z_score_list,cmap=plt.get_cmap('Reds'),marker='s')


    result_df_2=hit_df_sub.loc[(hit_df_sub['scaled_binding_score']>0.5)&((hit_df_sub['scaled_disrupt_score']>0.7)|(hit_df_sub['scaled_disrupt_score']<-0.7))]
    for index in result_df_2.index:

        a=result_df_2.loc[index]['scaled_binding_score']
        b=result_df_2.loc[index]['scaled_disrupt_score']
        c=result_df_2.loc[index]['motif'].split("__")[1]
        plt.text(float(b)+0.05, float(a)+0.00, c, fontsize=12)
    plt.plot([-1.3,1.3],[0.5,0.5],'k--',lw=1)
    plt.ylabel('Binding score')
    plt.xlabel('Disruption score')
    plt.xlim(-1.3,1.3)
    plt.colorbar()
    plt.savefig(pdffilename)

def plot_motif_specific_main(snp_motif_result_file, bg_score_folder, motif_id_name, pngfilename):
    snp_motif_result_df=pd.read_csv(snp_motif_result_file, sep='\t')
    bg_filename=os.path.join(bg_score_folder,motif_id_name+".scores")
    SNP_score_dataframe=pd.read_csv(bg_filename,sep="\t")
    SNP_score_dataframe.set_index('UID', inplace=True)
    score_df=SNP_score_dataframe[(SNP_score_dataframe['binding_ref']>0)|(SNP_score_dataframe['binding_alt']>0)]
    list_b=score_df['scaled_binding_score'] 
    list_a=score_df['scaled_disrupt_score']
    hit_df_sub=snp_motif_result_df.loc[snp_motif_result_df['motif']==motif_id_name]
    bg_bind_scaled_score_list=hit_df_sub['scaled_binding_score']
    bg_disrupt_scaled_score_list=hit_df_sub['scaled_disrupt_score']

    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w', edgecolor='w')
    ax = fig.add_subplot(111)
    plt.grid()
    #plt.plot([0, 1], [0, 0], 'k')
    #plt.plot([0, 0], [-1, 1], 'k')
    ##plot the points and then at the end
    ax.set_aspect('equal')
    plt.scatter(list_a, list_b, s=50,color='black',alpha=0.004)
    plt.scatter(bg_disrupt_scaled_score_list, bg_bind_scaled_score_list,  s=50,color='red')

    plt.ylabel('Binding score')
    plt.xlabel('Disruption score')
    plt.savefig(pngfilename)


    #plot area score distribution
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    sns.distplot( score_df["scaled_area_score"] , color="black", bins=50, label="SNP-background", kde_kws={"bw":0.02})
    sns.distplot( hit_df_sub["scaled_area_score"] , color="red", bins=50, label="SNP-hits", kde_kws={"bw":0.02})
    sns.plt.legend()
    sns.plt.xlim(-1,1)
    sns.plt.grid()
    axs.set(xlabel='Area Score')
    area_pdf_1=pngfilename.replace(".png","_area1.pdf")
    sns.plt.savefig(area_pdf_1)

    f, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    sns.distplot( score_df["scaled_area_score"] , color="black", bins=50, ax=axes[0],label="SNP-background")
    sns.distplot( hit_df_sub["scaled_area_score"] , color="red", bins=50, ax=axes[1], label="SNP-hits")
    sns.plt.xlim(-1,1)
    axes[0].legend()
    axes[0].set(xlabel='')
    axes[0].grid(True)
    axes[1].legend()
    axes[1].set(xlabel='Area Score')
    axes[1].grid(True)
    area_pdf_2=pngfilename.replace(".png","_area2.pdf")
    sns.plt.savefig(area_pdf_2)



def plot_motif_scattering_main(motiffile,all_expression_file,output_motif_summary_filename,motiffigure_filename,exp_cutoff_min, exp_cutoff_max, pvalue_cutoff_min,pvalue_cutoff_max):

    result_pvalue=pd.read_csv(motiffile,sep="\t")
    all_expression=pd.read_csv(all_expression_file,sep="\t")

    expression_zscore=[]
    expression_pvalue=[]
    count=0
    for index in result_pvalue.index:
        expression=result_pvalue.loc[index]['expression']
        motif_id=result_pvalue.loc[index]['motif']
        expression_list=all_expression.loc[all_expression['motif']==motif_id]
        expression_value_list0=expression_list[expression_list.columns[1:expression_list.shape[1]]].values.tolist()[0]
        expression_value_list = [x for x in expression_value_list0 if str(x) != 'nan']
        m=np.mean(expression_value_list)
        s=np.std(expression_value_list)
        if s==0:
            s=1
        z_score=(expression-m)/s
        p_value=scipy.stats.norm.cdf(z_score)
        expression_zscore.append(z_score)
        expression_pvalue.append(p_value)
        count=count+1

    result_pvalue['expression_zscore']=expression_zscore
    result_pvalue['expression_pvalue']=expression_pvalue

    sub_result_pvalue=result_pvalue.loc[((result_pvalue['expression_pvalue'] >= exp_cutoff_min) & (result_pvalue['expression_pvalue'] <= exp_cutoff_max))].sort_values(by=['expression_pvalue'],ascending=False)
    sub_result_pvalue=sub_result_pvalue.sort_values(by=['disrupt_pvalue'],ascending=True)
    pvalues=sub_result_pvalue['disrupt_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    sub_result_pvalue['tested_disrupt_FDR']=qvalue

    sub_result_pvalue=sub_result_pvalue.sort_values(by=['enhance_pvalue'],ascending=True)
    pvalues=sub_result_pvalue['enhance_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    sub_result_pvalue['tested_enhance_FDR']=qvalue

    sub_result_pvalue=sub_result_pvalue.sort_values(by=['abs_area_score_pvalue'],ascending=True)
    pvalues=sub_result_pvalue['abs_area_score_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    sub_result_pvalue['tested_abs_area_score_FDR']=qvalue

    dfright=sub_result_pvalue[['motif','tested_disrupt_FDR','tested_enhance_FDR','tested_abs_area_score_FDR']]


    result_pvalue=result_pvalue.sort_values(by=['disrupt_pvalue'],ascending=True)
    pvalues=result_pvalue['disrupt_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    result_pvalue['disrupt_FDR']=qvalue

    result_pvalue=result_pvalue.sort_values(by=['enhance_pvalue'],ascending=True)
    pvalues=result_pvalue['enhance_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    result_pvalue['enhance_FDR']=qvalue

    result_pvalue=result_pvalue.sort_values(by=['abs_area_score_pvalue'],ascending=True)
    pvalues=result_pvalue['abs_area_score_pvalue']
    qvalue=multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=True)[1]
    result_pvalue['abs_area_score_FDR']=qvalue

    merged_result_pvalue=pd.merge(result_pvalue, dfright, how='left', on='motif')
    merged_result_pvalue['minFDR']=merged_result_pvalue[['disrupt_FDR','enhance_FDR','abs_area_score_FDR','tested_disrupt_FDR','tested_enhance_FDR','tested_abs_area_score_FDR']].min(axis=1)
    merged_result_pvalue['minFDR_column']=merged_result_pvalue[['disrupt_FDR','enhance_FDR','abs_area_score_FDR','tested_disrupt_FDR','tested_enhance_FDR','tested_abs_area_score_FDR']].idxmin(axis=1)
    merged_result_pvalue=merged_result_pvalue.sort_values(by=['minFDR'],ascending=True)

    merged_result_pvalue.to_csv(output_motif_summary_filename,sep="\t",index=None)
    


    #finish preparing the whole motif file
    result_pvalue=pd.read_csv(output_motif_summary_filename,sep="\t",header=0)

    #sub_result_pvalue=result_pvalue.loc[((result_pvalue['expression_pvalue'] >= exp_cutoff_min) & (result_pvalue['expression_pvalue'] <= exp_cutoff_max))].sort_values(by=['expression_pvalue'],ascending=False)

    plt.figure(figsize=(8,8))
    x=-np.log10(result_pvalue['minFDR']+1E-50)
    y=result_pvalue['expression_pvalue']
    colors=[]
    for index in result_pvalue.index:
        colors.append('grey')
    plt.scatter(x, y, c=colors, alpha=0.3, s=50)


    sub_result_pvalue=result_pvalue.loc[((result_pvalue['expression_pvalue'] >= exp_cutoff_min) & (result_pvalue['expression_pvalue'] <= exp_cutoff_max))].sort_values(by=['expression_pvalue'],ascending=False)
    result_df_1=sub_result_pvalue.loc[(sub_result_pvalue['minFDR']<pvalue_cutoff_max)&(sub_result_pvalue['minFDR'] >= pvalue_cutoff_min)]
    x1=-np.log10(result_df_1['minFDR']+1E-50)
    y1=result_df_1['expression_pvalue']

    color1map={"disrupt_FDR":"blue","enhance_FDR":"red","abs_area_score_FDR":"orange","tested_disrupt_FDR":"blue","tested_enhance_FDR":"red","tested_abs_area_score_FDR":"orange"}
    colors1=[color1map[key] for key in result_df_1['minFDR_column']]

    '''
    for index in result_df_1.index:
        if result_df_1.loc[index]['leftside_pvalue']<result_df_1.loc[index]['rightside_pvalue']:
            colors1.append('red')
        else:
            colors1.append('blue')
        if result_df_1.loc[index]['correct_pvalue']==0:
            continue
        a=-np.log10(result_df_1.loc[index]['correct_pvalue'])
        b=result_df_1.loc[index]['expression_pvalue']
        c=result_df_1.loc[index]['motif'].split("__")[1]
        #plt.text(float(a)+0.1, float(b), c, fontsize=8)
    '''
    plt.scatter(x1, y1, c=colors1, alpha=0.9, s=50)

    plt.xlabel('-log10 corrected p-value')
    plt.ylabel('Specificity')
    plt.title('Expression vs enriched p-value for Transcription Factors')
    plt.savefig(motiffigure_filename)


    zoomin_figure_filename=motiffigure_filename.replace(".pdf","")+"_zoom.pdf"
    plt.figure(figsize=(8,8))
    x1=-np.log10(result_df_1['minFDR']+1E-50)
    y1=result_df_1['expression_pvalue']
    
    colors1=[color1map[key] for key in result_df_1['minFDR_column']]
    
    for index in result_df_1.index:
        a=-np.log10(result_df_1.loc[index]['minFDR']+1E-50)
        b=result_df_1.loc[index]['expression_pvalue']
        c=result_df_1.loc[index]['motif'].split("__")[1]
        plt.text(float(a)-0.2, float(b)+0.005, c, fontsize=16)
    
    plt.scatter(x1, y1, c=colors1, alpha=0.9, s=50)

    plt.xlabel('-log10 corrected p-value')
    plt.ylabel('Specificity')
    plt.title('Expression vs enriched p-value for Transcription Factors')

    plt.savefig(zoomin_figure_filename)
