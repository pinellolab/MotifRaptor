# MotifRaptor

## motivation

The current challenges in post-GWAS study is to infer the potential mechanism and filter further for potential causal SNPs. The current methods are either requiring too much extra experimental data or lacking efficient and systematic way to do it. We would like to design a toolkit that simply an universal sequence-based motif scanning for trait associated SNPs that has been suffering from computational complexity and false positives for years.

## Design

Motif-raptor is a toolkit that solely depending on genome-wide motif scanning and thus can be universally applied in any phenotype related functional genomics. We designed a prefix-suffix based method to achieve the computational complexity of the genome wide screening. We also designed informative quantification scores and effective statistics to characterize important cell or tissue types, TFs, and SNP sites given GWAS summary statistics. We tested our method on Rheumatoid Arthritis and re-discover TFs and SNPs specifically contributing to the trait in immune cells. We provide several demos on this site.


## Installation

1. Download files or git clone from this repository and keep them in the same folder

2. Extra python packages:
   ```
    pip install twobitreader
    pip install pybedtools
   ```
3. Install mksary

   3.0 Check the file "mksary" in the folder MotifRaptor/SNPScanner, if it's runnable, you don't need to build your own.

   3.1 Download [libdivsufsort](https://goo.gl/hUjvMF), which is the library available [in GitHub](https://github.com/y-256/libdivsufsort) changed in the file "mksary.c" which now includes the computation of the LCP array.
   
   3.2 Build   
   ```
    cd libdivsufsort 
    mkdir build 
    cd build
    make 
    sudo make install
   ```

    3.3 Copy the executable "mksary" into your python-working directory, here it should be MotifRaptor/SNPScanner.

4. Compile the cython code: cythonize -a -i motif_matching_lcp.pyx

## Database

1. Download the database from Dropbox link. This database contains essential data for general analysis, including DHS tracks, TF RNA-seq expressions, TF motifs, and TF pre-calucated scores. 

2. Unzip the Database.zip, and this folder should be inside the main folder of the MotifRaptor package.

## Example

### Data

The data are in folder Demo/RA. There is a jupyter notebook can directly guide the steps.

### step 0. pre-processing
This step is to generate standard format for hit SNPs and background SNPs. We need three standard format from GWAS summary statistics. Based on the p-value cutoff, the user should define the SNP hits and SNP non-hits lists in text files, and a VCF file for SNP hits. These file should be already in the folder Demo/RA.

   ```
    cd Demo/RA 
    ll 
    hitSNP_list.txt
    nonhitSNP_list.txt 
    hitSNP_list.vcf
   ```
If you don't see these files, you can always make your own with a short python script. This needs to be customized base on each different format in GWAS summary statistics, and applying different cut-offs as you like.
   ```
    import numpy as np
    import pandas as pd
    sumstatsfile="RA_GWASmeta_TransEthnic_v2.txt"
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


   ```
### step1. run cell type or tissue type characterization

   ```
    python package_path/MotifRaptor/MotifRaptor.py celltype -vcf hitSNP_list.vcf \
        -sh hitSNP_list.txt -sn nonhitSNP_list.txt -wd step1_out/ -p 2
    
   ```
   
   Optional trial:
   
   Motif raptor also provide module interface to run it in your own python code or Jupyter notebook for lightweight task. For example, you can re-plot the figure after running the previous steps.
   
   ```
    #optional: plot figure using module CellTypeAnlayzer after running the command
    from CellTypeAnalyzer import CellTypeAnalysis
    CellTypeAnalysis.plotfigure_main("all_sorted.pvalue","plot_all_cell_type.pdf")
   ```


<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic1.png" alt="drawing" width="800"/>

### step2. run motif discovery and filtering

   ```
    #run through the motif scan on all the SNPs
     python package_path/MotifRaptor/MotifRaptor.py snpmotif \
     -wd step2_out \
     -c "CD8-positive, alpha-beta T cell" \
     -sb step1_out/hitSNP_DHSunion_list.bed \
     -se step1_out/hit_snps_seq_pickle.df.txt \
     -bg genome -m test_motif.txt -p 2
   ```
    
   ```
    #run through the motif scan on all the SNPs
      python package_path/MotifRaptor/MotifRaptor.py motiffilter \
      -wd step2_out/motif_result \
      -ms step2_out/motif_result/all_motifs.pvalue
   ```
   
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic2.png" alt="drawing" width="800"/>

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic3.png" alt="drawing" width="800"/>


   
### step3. run filter for SNPs and motifs and draw plots
    
   ```
    # per motif analysis
     python package_path/MotifRaptor/MotifRaptor.py motifspecific \
     -wd step3_out \
     -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
     -md MA0105.1__NFKB1 \
     -bs step2_out/motif_result/background_files/ 
   
   ```
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic6.png" alt="drawing" width="800"/>

   ```
    # per SNP analysis
    python package_path/MotifRaptor/MotifRaptor.py snpspecific \
    -wd step3_out \
    -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
    -snp rs7528684 
   
   ```
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic4.png" alt="drawing" width="800"/>
    
   ```
    # draw radar plot for instereing motif and SNP events
    python package_path/MotifRaptor/MotifRaptor.py snpmotifradar \
    -wd step3_out \
    -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
    -pid rs7528684:MA0105.1__NFKB1
   ```

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic5.png" alt="drawing" width="800"/>



## Modules and usages

   ```
    usage: MotifRaptor [-h] [--version]
                   {celltype,snpmotif,motiffilter,motifspecific,snpspecific,snpmotifradar}

Analyze motifs and SNPs in the dataset.

positional arguments:
  {celltype,snpmotif,motiffilter,motifspecific,snpspecific,snpmotifradar}
                        help for subcommand: celltype, snpmotif, motiffilter,
                        motifspecific, snpspecific
    celltype            cell type or tissue type analysis help
    snpmotif            snp motif test help
    motiffilter         motifs filtering help
    motifspecific       motifs specific analysis help
    snpspecific         SNP specific analysis help
    snpmotifradar       SNP motif radar plot help

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
    
   ```
    
## Optional: Building disruption scores genome wide for a list of SNPs from a VCF file

### Data
1. VCF file:

   Download VCF file from [here](https://www.dropbox.com/s/9gztf4mdblc44jo/1000G.EUR.QC.plink.simple.vcf?dl=0)

2. Genome 2Bit File:
   
   Download from UCSC website [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) 

### Steps

1. Extract the sequences starting from a reference genome and a VCF file:
     ```
      python Process_vcf.py
     ```

2. Index the obtaines sequences with the variants:
    
     ```
      python Motif_Scan.py index -g genome_database -p number_of_threads
     ```
    
3. Scan motifs using PFM files and the obtained index: 

     ```
      python Motif_Scan.py pfmscan -gi indexed_genome_database -pfm motif_pfm_folder -mo motifscan_result -p number_of_threads
     ```

