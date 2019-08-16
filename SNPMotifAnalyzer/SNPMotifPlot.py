import numpy as np
import numbers
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

def qqplot(x, y, quantiles=None, interpolation='nearest', ax=None, rug=False,
           rug_length=0.05, rug_kwargs=None, **kwargs):
    """Draw a quantile-quantile plot for `x` versus `y`.

    Parameters
    ----------
    x, y : array-like
        One-dimensional numeric arrays.

    ax : matplotlib.axes.Axes, optional
        Axes on which to plot. If not provided, the current axes will be used.

    quantiles : int or array-like, optional
        Quantiles to include in the plot. This can be an array of quantiles, in
        which case only the specified quantiles of `x` and `y` will be plotted.
        If this is an int `n`, then the quantiles will be `n` evenly spaced
        points between 0 and 1. If this is None, then `min(len(x), len(y))`
        evenly spaced quantiles between 0 and 1 will be computed.

    interpolation : {‘linear’, ‘lower’, ‘higher’, ‘midpoint’, ‘nearest’}
        Specify the interpolation method used to find quantiles when `quantiles`
        is an int or None. See the documentation for numpy.quantile().

    rug : bool, optional
        If True, draw a rug plot representing both samples on the horizontal and
        vertical axes. If False, no rug plot is drawn.

    rug_length : float in [0, 1], optional
        Specifies the length of the rug plot lines as a fraction of the total
        vertical or horizontal length.

    rug_kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.axvline() and
        matplotlib.axes.Axes.axhline() when drawing rug plots.

    kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.scatter() when drawing
        the q-q plot.
    """
    # Get current axes if none are provided
    if ax is None:
        ax = plt.gca()

    if quantiles is None:
        quantiles = min(len(x), len(y))

    # Compute quantiles of the two samples
    if isinstance(quantiles, numbers.Integral):
        quantiles = np.linspace(start=0, stop=1, num=int(quantiles))
    else:
        quantiles = np.atleast_1d(np.sort(quantiles))
    #x_quantiles = np.quantile(x, quantiles, interpolation=interpolation)
    a=pd.Series(x.tolist())
    x_quantiles = a.quantile(quantiles, interpolation=interpolation)
    #y_quantiles = np.quantile(y, quantiles, interpolation=interpolation)
    b=pd.Series(y.tolist())
    y_quantiles = b.quantile(quantiles, interpolation=interpolation)
    # Draw the rug plots if requested
    if rug:
        # Default rug plot settings
        rug_x_params = dict(ymin=0, ymax=rug_length, c='gray', alpha=0.5)
        rug_y_params = dict(xmin=0, xmax=rug_length, c='gray', alpha=0.5)

        # Override default setting by any user-specified settings
        if rug_kwargs is not None:
            rug_x_params.update(rug_kwargs)
            rug_y_params.update(rug_kwargs)

        # Draw the rug plots
        for point in x:
            ax.axvline(point, **rug_x_params)
        for point in y:
            ax.axhline(point, **rug_y_params)

    # Draw the q-q plot
    ax.scatter(x_quantiles, y_quantiles, **kwargs)

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

def plot_motif_snp_pair_main(snp_motif_result_file, all_expression_file,snp_id, motif_id_name, rsid_motifid, pdffilename, snpfeaturefile, catofile):
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

    #old way to plot cato score
    cato_df=pd.read_csv(catofile,sep="\t",header=0)
    motif_id=motif_id_name.split("__")[0]
    cato_score_hits=cato_df[(cato_df['rsid']==snp_id)&(cato_df['motifid']==motif_id)]['score'].values
    cato_score=0.0
    if len(cato_score_hits)>0:
        cato_score=cato_score_hits[0]
    radar_dict['CATO']=cato_score

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
    #sns.plt.legend()
    #sns.plt.xlim(-1,1)
    #sns.plt.grid()
    plt.legend()
    plt.xlim(-1,1)
    plt.grid()
    axs.set(xlabel='Area Score')
    area_pdf_1=pngfilename.replace(".png","_area1.pdf")
    #sns.plt.savefig(area_pdf_1)
    plt.savefig(area_pdf_1)

    f, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    sns.distplot( score_df["scaled_area_score"] , color="black", bins=50, ax=axes[0],label="SNP-background")
    sns.distplot( hit_df_sub["scaled_area_score"] , color="red", bins=50, ax=axes[1], label="SNP-hits")
    #sns.plt.xlim(-1,1)
    plt.xlim(-1,1)
    axes[0].legend()
    axes[0].set(xlabel='')
    axes[0].grid(True)
    axes[1].legend()
    axes[1].set(xlabel='Area Score')
    axes[1].grid(True)
    area_pdf_2=pngfilename.replace(".png","_area2.pdf")
    #sns.plt.savefig(area_pdf_2)
    plt.savefig(area_pdf_2)

    #plot Q-Q plot for two samples
    X_array=np.array(score_df["scaled_area_score"])
    Y_array=np.array(hit_df_sub["scaled_area_score"])
    qqplot_pdf=pngfilename.replace(".png","_qqplot.pdf")
    fig=plt.figure(figsize=(10, 10), dpi=80)
    ax = fig.add_subplot(111)
    #a=X_array[1:2000]
    #b=X_array[4001:6000]
    #qqplot(a,b,c='r', alpha=0.5, edgecolor='k')
    qqplot(X_array, Y_array, c='r', alpha=0.5, edgecolor='k',quantiles=20)
    plt.plot([-1,1],[-1,1],'b--')
    plt.xlabel('SNP-background Area Score')
    plt.ylabel('SNP-hits Area Score')
    plt.tight_layout()
    #plt.title('Q-Q plot')
    plt.savefig(qqplot_pdf)

    X_all=score_df["scaled_area_score"].tolist()
    ind=np.random.choice(len(X_all), min(20000,len(X_all)), replace=False).tolist()
    X_list=[X_all[index] for index in ind]
    Y_list=hit_df_sub["scaled_area_score"].tolist()
    #ind=np.random.choice(len(Y_all), min(20000,len(Y_all)), replace=False).tolist()
    #Y_list = [Y_all[index] for index in ind]
    X_list_len=len(X_list)
    Y_list_len=len(Y_list)
    SNP_list=['background']*X_list_len+['hits']*Y_list_len
    Motif_list=[motif_id_name]*(X_list_len+Y_list_len)
    df_violinplot=pd.DataFrame([],index=range(0,X_list_len+Y_list_len), columns=['Area Score','Motif','SNP'])
    df_violinplot['Area Score']=X_list+Y_list
    df_violinplot['Motif']=Motif_list
    df_violinplot['SNP']=SNP_list
    violinplot_pdf=pngfilename.replace(".png","_violinplot.pdf")
    fig=plt.figure(figsize=(10, 10), dpi=80)
    ax = fig.add_subplot(111)
    palette ={"hits":"red", "background":"grey"}
    ax = sns.violinplot(x="Motif", y="Area Score", hue="SNP", data=df_violinplot, palette=palette, split=True, scale="area", inner="quartile", scale_hue=True,linewidth=1, bw=0.2)
    plt.tight_layout()    
    plt.savefig(violinplot_pdf)


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
