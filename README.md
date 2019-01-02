# MotifRaptor

## motivation

The current challenges in post-GWAS study is to infer potential mechanism and narrow down the validation list for potential causal SNPs. The current methods are either requiring too much extra experimental data or lacking efficient and systematic way to do it. We would like to design a toolkit performing an universal sequence-based motif scanning for trait associated SNPs that has been suffering from computational complexity and false positives for years.

## Design

Motif-Raptor (**Motif** dis**R**uption **A**ssociated **P**olymorphism for **T**ranscription fact**OR**s) is a toolkit that solely depending on genome-wide motif scanning around SNPs and thus can be universally applied in any phenotype related functional genomics. We designed a prefix-suffix based method to conqure the computational complexity of this genome wide screening. We also designed informative quantification scores and effective statistics to characterize important cell or tissue types, prioritize Transcription Factors and SNP sites given GWAS summary statistics. We tested our method on Rheumatoid Arthritis and re-discover TFs and SNPs specifically contributing to the trait in immune cells.



## Installation

1. Download files or git clone from this repository and keep them in the same folder
   ```
   git clone https://github.com/pinellolab/MotifRaptor.git
   tar zvxf MotifRaptor.tar.gz
   ```

2. Extra python packages:
   ```
    pip install twobitreader
    pip install pybedtools
   ```
3. Install mksary

   3.0 Check the file "mksary" in the folder MotifRaptor/SNPScanner, if it's runnable, you don't need to build your own and go to step 4.
   ```
   cd MotifRaptor/SNPScanner
   ./mksary --help
   ```

   3.1 Download [libdivsufsort](https://goo.gl/hUjvMF), which is the library available [in GitHub](https://github.com/y-256/libdivsufsort) changed in the file "mksary.c" which now includes the computation of the LCP array.
   
   3.2 Build   
   ```
    cd libdivsufsort 
    mkdir build 
    cd build
    make 
    sudo make install
   ```

    3.3 Copy the executable "mksary" into your python-working directory, here it should be MotifRaptor/SNPScanner. Double check if it's working now.
   ```
   ./mksary --help
   ```    

4. Compile the cython code: cythonize -a -i motif_matching_lcp.pyx

## Database

1. Download the database from Dropbox link. This database contains essential data for general analysis, including DHS tracks, TF RNA-seq expressions, TF motifs, and TF pre-calucated scores. 
   ```
    wget ****
   ```

2. Unzip the Database.zip, and this folder should be inside the main folder of the MotifRaptor package.
   ```
   mv Database.zip MotifRaptor/
   cd MotifRaptor
   unzip Database.zip
   ```
## Example

### Data

The data are in folder Demo/RA. There is a jupyter notebook can directly guide you these steps.

### step 0. pre-processing
This step is to generate standard format for hit SNPs and background SNPs. We need three standard format from GWAS summary statistics. Based on the p-value cutoff, the user should define the SNP hits and SNP non-hits lists in text files, and a VCF file for SNP hits. These file should be already in the folder Demo/RA.

   ```
    cd Demo/RA 
    ll 
    hitSNP_list.txt
    nonhitSNP_list.txt 
    hitSNP_list.vcf
   ```
If you don't see these files, you can always make your own with a short python script. This needs to be customized base on each different format in GWAS summary statistics (Okada et al. 2010 Nature), and applying different cut-offs as you like.
   ```
   wget https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Okada/RA_GWASmeta_TransEthnic_v2.txt.gz
   gunzip RA_GWASmeta_TransEthnic_v2.txt.gz
   
   python
   
   >
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
   >

   ```
The outcome of these codes are generate SNP list for GWAS hit and non-hits based on some cut-off. 

   
### step1. run cell type or tissue type characterization

   ```
    python package_path/MotifRaptor/MotifRaptor.py celltype -vcf hitSNP_list.vcf \
        -sh hitSNP_list.txt -sn nonhitSNP_list.txt -wd step1_out/ -p 2
    
   ```
**Input:** *hitSNP_list.txt* and *nonhitSNP_list.txt* are two files with the following format:

ID | CHR | POS
------------ | ------------- | ------------- 
rs2258734 | 1 | 2483961

*hitSNP_list.vcf* is the file with the following format (with two more columns of the polymorphism information:

ID | CHR | POS | REF | ALT
------------ | ------------- | ------------- | ------------- | ------------- 
rs2258734 | 1 | 2483961 | A | G


   *Optional trial:*
   
   *Motif raptor also provide module interface to run it in your own python code or Jupyter notebook for lightweight task. For example, you can re-plot the figure after running the previous steps.*
   
   ```
    #optional: plot figure using module CellTypeAnlayzer after running the command
    from CellTypeAnalyzer import CellTypeAnalysis
    CellTypeAnalysis.plotfigure_main("all_sorted.pvalue","plot_all_cell_type.pdf")
   ```

**Output:** This visualization ranks the associated cell or tissue types by both p-values (bar length) and odd ratio (numbers behind the bar).

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic1.png" alt="drawing" width="600"/>

**Usage:**
   ```
   python MotifRaptor.py celltype --help
   
   usage: MotifRaptor celltype [-h] [-vcf SNP_HIT_VCF] [-sh SNP_HIT]
                            [-snh SNP_NON_HIT] [-wd WORKDIR] [-p THREAD_NUM]

optional arguments:
  -h, --help            show this help message and exit
  -vcf SNP_HIT_VCF, --snp_hit_withseq SNP_HIT_VCF
                        need header and columns in this text file with
                        sequence (CHR is only a number): ID CHR POS REF ALT
  -sh SNP_HIT, --snp_hit SNP_HIT
                        need header and columns in this text file (CHR is only
                        a number): ID CHR POS
  -snh SNP_NON_HIT, --snp_non_hit SNP_NON_HIT
                        need header and columns in this text file (CHR is only
                        a number): ID CHR POS
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -p THREAD_NUM, --threads THREAD_NUM
                        number of threads
   ```

### step2. run motif discovery and filtering

#### step2.1 run through statistics and essential analysis

   ```
    #run through the motif scan on all the SNPs
     python package_path/MotifRaptor/MotifRaptor.py snpmotif \
     -wd step2_out \
     -c "CD8-positive, alpha-beta T cell" \
     -sb step1_out/hitSNP_DHSunion_list.bed \
     -se step1_out/hit_snps_seq_pickle.df.txt \
     -bg genome -m test_motif.txt -p 2
   ```
    

**Input:** The bed file for SNP hits and the sequence information are obtained from Step1 "cell type characterization". But you need to determine the background and the motifs you want to test.

You may use "genome" in "-g" to use genome wide SNPs as the baseline distribution.
You may use "all" in "-m" to test all of the Transcription Factors collected in the database. But you may specificy a file *test_motif.txt* to run a test for only a few Transcription Factors. Each motif should take a line, with the format of "motifID\_\_motifname" which pwm files can be found in the Database. For example:
   ```
    MA0062.1__GABPA
    MA0105.1__NFKB1
    MA0518.1__Stat4
   ```
**Output:** Scores for SNP hits are calculated as text files, including a huge background folder. The outcome of this partial step is not viewable. You need to run the next step to get a summazied text and figures.

   ```
   python MotifRaptor.py snpmotif --help
usage: MotifRaptor snpmotif [-h] [-wd WORKDIR] [-c CELL_TYPE]
                            [-sb HIT_SNP_UNION_BED] [-se HIT_SNP_UNION]
                            [-bg BG_SNPS] [-m MOTIF_LIST] [-p THREAD_NUM]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -c CELL_TYPE, --cell_type CELL_TYPE
                        Cell type or Tissue type ID
  -sb HIT_SNP_UNION_BED, --snp_hit_bed HIT_SNP_UNION_BED
                        hit snps on union bed file, usually from step1
  -se HIT_SNP_UNION, --snp_hit_seq HIT_SNP_UNION
                        hit snps with sequence information, usually from step1
  -bg BG_SNPS, --bg_snp BG_SNPS
                        background snp list file or (genome)
  -m MOTIF_LIST, --motifs MOTIF_LIST
                        motif list file, no header, or (all)
  -p THREAD_NUM, --threads THREAD_NUM
                        number of threads
   ```

#### step2.2 apply the summary, filter and plot figures
   ```
    #run through the motif scan on all the SNPs
      python package_path/MotifRaptor/MotifRaptor.py motiffilter \
      -wd step2_out/motif_result \
      -ms step2_out/motif_result/all_motifs.pvalue
   ```
   
**Input:** Just specify the working directory from previous step.
   
**Output:** The visualization shows the distribution and scoring for the Transcription Factors and narrow down (zoomed in) with the significant/interesting ones.

*all motifs*

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic2.png" alt="drawing" width="400"/>

*zoomed in motifs*

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic3.png" alt="drawing" width="400"/>

   ```
   python MotifRaptor.py motiffilter --help
usage: MotifRaptor motiffilter [-h] [-wd WORKDIR] [-ms MOTIFFILE]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -ms MOTIFFILE, --motif_summary MOTIFFILE
                        Motif Summary File, usually from step2
   ```
   
### step3. run filter for SNPs and motifs and draw plots

#### step3.1 TF specific plot
   ```
    # per motif analysis
     python package_path/MotifRaptor/MotifRaptor.py motifspecific \
     -wd step3_out \
     -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
     -md MA0105.1__NFKB1 \
     -bs step2_out/motif_result/background_files/ 
   
   ```
**Input:** Just specify the output directory from step2.
   
**Output:** Scatter plot for both SNP hits and non-hits for a picked motif.
   
   
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic6.png" alt="drawing" width="400"/>

#### step3.2 SNP specific plot
   ```
    # per SNP analysis
    python package_path/MotifRaptor/MotifRaptor.py snpspecific \
    -wd step3_out \
    -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
    -snp rs7528684 
   
   ```
   
**Input:** Just specify the output directory from step2.
   
**Output:** Scatter plot for transcription factors for a picked SNP.
  
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic4.png" alt="drawing" width="400"/>

#### step3.3 SNP-TF plot    
   ```
    # draw radar plot for instereing motif and SNP events
    python package_path/MotifRaptor/MotifRaptor.py snpmotifradar \
    -wd step3_out \
    -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
    -pid rs7528684:MA0105.1__NFKB1
   ```

**Input:** Just specify the output directory from step2.
   
**Output:** Radar plot for the features/scores for interesting SNP-motif events.


<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic5.png" alt="drawing" width="400"/>



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


