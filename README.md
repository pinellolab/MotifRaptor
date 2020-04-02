# MotifRaptor

## Overview

##### Motivation: 
Genome-wide association studies (GWAS) have identified thousands of common trait-associated genetic variants but interpretation of their function remains challenging. These genetic variants can overlap the binding sites of transcription factors (TFs) and therefore could alter gene expres-sion. However, we currently lack a systematic understanding on how this mechanism contrib-utes to phenotype. 
##### Results: 
We present Motif-Raptor, a TF-centric computational tool that integrates sequence-based predic-tive models, chromatin accessibility, gene expression datasets and GWAS summary statistics to systematically investigate how TF function is affected by genetic variants. Given trait associated non-coding variants, Motif-Raptor can recover relevant cell types and critical TFs to drive hy-potheses regarding their mechanism of action. We tested Motif-Raptor on complex traits such as rheumatoid arthritis and red blood cell count and demonstrated its ability to prioritize relevant cell types, potential regulatory TFs and non-coding SNPs which have been previously characterized and validated.



## Installation

1. Install from Bioconda

   ```
   conda create -n motifraptor_env python=3.6
   source activate motifraptor_env
   conda install -c bioconda motifraptor
   ```
   
2. Simple test

   Activate the conda environment before running the program.
   ```
   source activate motifraptor_env
   ```
   Motif Raptor supports two different ways to run it.
   
   (1) Run from the command line (recommended)
   ```
   MotifRaptor --version
   ```
   (2) Load as a module
   ```
   python
   >>>import MotifRaptor
   >>>MotifRaptor.__version__
   ```
   If you see the version number, congratulations!
   
## Motif-Raptor Modules Overview

   ```
   MotifRaptor --help
   ```
   ```
usage: MotifRaptor [-h] [--version]
                   
                   {preprocess,preprocess_ukbb_v3,celltype,snpmotif,snpfeature,motiffilter,motifspecific,snpspecific,snpmotifradar,snpindex,snpscan,set,info}
                   ...

Analyze motifs and SNPs in the dataset.

positional arguments:
  {preprocess,preprocess_ukbb_v3,celltype,snpmotif,snpfeature,motiffilter,motifspecific,snpspecific,snpmotifradar,snpindex,snpscan,set,info}
                        help for subcommand: celltype, snpmotif, snpfeature,
                        motiffilter, motifspecific, snpspecific
    preprocess          Pre-process the summary statistics
    preprocess_ukbb_v3  Pre-process the summary statistics from UKBB version 3
                        TSV files
    celltype            cell type or tissue type analysis help
    snpmotif            snp motif test help
    snpfeature          snp feature help
    motiffilter         motifs filtering help
    motifspecific       motifs specific analysis help
    snpspecific         SNP specific analysis help
    snpmotifradar       SNP motif radar plot help
    snpindex            index the SNPs (with flanking sequences) help
    snpscan             scan SNP database (already indexed) help
    set                 Set Path and Global Values
    info                Get Informationa and Print Global Values

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
    
   ```   
<!--   
2. Download and Install Database

Download the database from Dropbox link. This database contains essential data for general analysis, including DHS tracks, TF RNA-seq expressions, TF motifs, and TF pre-calucated scores. 

*The following file is a database with a small number of TFs, in order to test the tutorial example on a regular machine.*

   ```
   wget https://www.dropbox.com/s/gxeyzgl5m0u55w8/Database.zip
   unzip Database.zip
   ```
Move database into the package folder
   ```
   #If a different python version was installed, change "python3.6" into a correspondent name

   mv Database $CONDA_PREFIX/lib/python3.6/site-packages/MotifRaptor/
   
   #or instead, you can generate a link to link from the package installation path to the database to avoid moving files
   
   cd $CONDA_PREFIX/lib/python3.6/site-packages/MotifRaptor/
   ln -s pathTo/Database Database

   ```
 *For a complete TF list, please download here. We strongly recommend run the program on a server. *
 
   ```
   wget https://www.dropbox.com
   unzip Database.zip
   ```
-->


## Tutorial


### step 0. prepare input data (pre-processing)

Run Motif-Raptor from GWAS summary statistics. You may get summary statistics from UKBiobank, published paper, or other resources. These files may provide diffrent information. Please make sure the file contains the following columns.

ID | CHR | POS | REF | ALT | p-value
------------ | ------------- | ------------- | ------------- | -------------  | ------------- 
rs2258734 | 1 | 2483961 | A | G | 0.003

From these columns and the following code, you can make three files *SNP hits* and *SNP non-hits* lists in text files, and *a VCF file for SNP hits*

```
python code
```


For a simple example, please download the following three files directly. These files are generated from (Okada et al. 2010 Nature)

```

wget https://www.dropbox.com/s/as5i16fur4g90m0/hitSNP_list.txt
wget https://www.dropbox.com/s/pq4ln4x19k567ua/nonhitSNP_list.txt
wget https://www.dropbox.com/s/939b3eb8cu1pd56/hitSNP_list.vcf

```
**Check Input Format:** *hitSNP_list.txt* and *nonhitSNP_list.txt* are two files with the following format:

ID | CHR | POS
------------ | ------------- | ------------- 
rs2258734 | 1 | 2483961

*hitSNP_list.vcf* is the file with the following format (with two more columns of the polymorphism information:

ID | CHR | POS | REF | ALT
------------ | ------------- | ------------- | ------------- | ------------- 
rs2258734 | 1 | 2483961 | A | G

<!--
 *Optional step:*
 
* Although GWAS summary statistics files may have different columns, given a data file with at least the following columns you can always make the above three input files on your own, with a simple python script.

ID | CHR | POS | REF | ALT | p-value
------------ | ------------- | ------------- | ------------- | -------------  | ------------- 
rs2258734 | 1 | 2483961 | A | G | 0.003


*For example, download the original data file from (Okada et al. 2010 Nature), and applying your own cut-offs to define hits and nonhits.*
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
The outcome of these codes are the three input files mentioned above. 
-->
   
### step1. run cell type or tissue type characterization


   
**Example:** 
   ```
   MotifRaptor celltype -vcf hitSNP_list.vcf -sh hitSNP_list.txt -sn nonhitSNP_list.txt -wd step1_out/ -p 3
   ```


**Input:** Three input files from step0.
**Output:** This visualization ranks the associated cell or tissue types by both p-values (bar length) and odd ratio (numbers behind the bar). Original text file is "*all_sorted.pvalue*"

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic1.png" alt="drawing" width="600"/>

   *Optional step:*
   
   *Motif raptor also provide module interface to run it in your own python code or Jupyter notebook for lightweight task. For example, you can re-plot the figure after running the previous steps.*
   
   ```
    #optional: plot figure using module CellTypeAnlayzer after running the command
    from MotifRaptor.CellTypeAnalyzer import CellTypeAnalysis
    CellTypeAnalysis.plotfigure_main("step1_out/testcelltype/all_sorted.pvalue","plot_all_cell_type.pdf")
   ```


**Usage:**
   ```
   MotifRaptor celltype --help
   
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
   MotifRaptor snpmotif -wd step2_out -c "CD8-positive, alpha-beta T cell" \
   -sb step1_out/hitSNP_DHSunion_list.bed -se step1_out/hit_snps_seq_pickle.df.txt -bg genome \
   -m all -p 3
   ```
    

**Input:** The bed file for SNP hits and the sequence information are obtained from Step1 "cell type characterization". But you need to determine the background and the motifs you want to test.

You may use "genome" in "-g" to use genome wide SNPs as the baseline distribution.
You may use "all" in "-m" to test all of the Transcription Factors collected in the database. But you may specificy a file *test_motif.txt* to run a test for only a few Transcription Factors. Each motif should take a line, with the format of "motifID\_\_motifname" which pwm files can be found in the Database. For example:
   ```
    MA0105.1__NFKB1
    MA0518.1__Stat4
   ```
**Output:** Scores for SNP hits are calculated as text files, including a huge background folder. The outcome of this partial step is not viewable. You need to run the next step to get a summazied text and figures.

**Usage:**
   ```
   MotifRaptor snpmotif --help
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

##### step2.2.1 apply filter on TF summary file and visualize the TF plots

**Example:**  
   
   ```
   #calculate statistics for each transcription factor and apply filtering (i.e. expression and pvalue)
   MotifRaptor motiffilter -wd step2_out/motif_result -ms step2_out/motif_result/all_motifs.pvalue
   ```
   
**Input:** Just specify the working directory from previous step.
   
**Output:** The visualization shows the distribution and scoring for the Transcription Factors and narrow down (zoomed in) with the significant/interesting ones.

*all motifs*

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic2.png" alt="drawing" width="400"/>

*zoomed in motifs*

<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic3.png" alt="drawing" width="400"/>

**Usage:**

   ```
   MotifRaptor snpfeature -h
   usage: MotifRaptor snpfeature [-h] [-wd WORKDIR] [-c CELL_TYPE]
                              [-cb SNP_BED_FILES]

 optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -c CELL_TYPE, --cell_type CELL_TYPE
                        Cell type or Tissue type ID
  -cb SNP_BED_FILES, --snp_bed_files SNP_BED_FILES
                        SNP cell type bed file folder, usually from step1
   ```
##### step2.2.2 pre-calculate the SNP-wise features to prepare visualization in the next steps
   ```
   #calculate SNP wise features (prepare for step3)
   MotifRaptor snpfeature -wd step2_out -c "CD8-positive, alpha-beta T cell" -cb step1_out/testcelltype/bedfiles
   ```
   
**Input:** Just specify the working directory from previous step.
   
**Output:** The SNP-wise annotations will be added to the folder.

**Usage:**

   ```
   MotifRaptor snpfeature -h
   usage: MotifRaptor snpfeature [-h] [-wd WORKDIR] [-c CELL_TYPE]
                              [-cb SNP_BED_FILES]

 optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -c CELL_TYPE, --cell_type CELL_TYPE
                        Cell type or Tissue type ID
  -cb SNP_BED_FILES, --snp_bed_files SNP_BED_FILES
                        SNP cell type bed file folder, usually from step1
   ```

### step3. SNP and TF events visualization

#### step3.1 TF specific plot
 **Example:**  
  ```
    # per motif analysis
   MotifRaptor motifspecific \
   -wd step3_out \
   -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
   -md MA0105.1__NFKB1 \
   -bs step2_out/motif_result/background_files/
   
   ```
**Input:** Just specify the output directory from step2.
   
**Output:** Scatter plot for both SNP hits and non-hits for a picked motif.
   
   
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic6.png" alt="drawing" width="400"/>

**Usage:**
```
usage: MotifRaptor motifspecific [-h] [-wd WORKDIR] [-sm SNP_MOTIF_FILE]
                                 [-md MOTIF_ID_NAME] [-bs BG_SCORE_FOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -sm SNP_MOTIF_FILE, --snp_motif_file SNP_MOTIF_FILE
                        SNP motif pair-wise list File, usually from step2
  -md MOTIF_ID_NAME, --motif_id MOTIF_ID_NAME
                        motif id with name, in the format of motifID__NAME
  -bs BG_SCORE_FOLDER, --bg_score_folder BG_SCORE_FOLDER
                        background score folder, usually from step2
```

#### step3.2 SNP specific plot
**Example:**  
   ```
    # per SNP analysis
   MotifRaptor snpspecific \
   -wd step3_out \
   -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
   -snp rs7528684 
   
   ```
   
**Input:** Just specify the output directory from step2.
   
**Output:** Scatter plot for transcription factors for a picked SNP.
  
<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic4.png" alt="drawing" width="400"/>

**Usage:**
```
usage: MotifRaptor snpspecific [-h] [-wd WORKDIR] [-sm SNP_MOTIF_FILE]
                               [-snp SNP_ID]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -sm SNP_MOTIF_FILE, --snp_motif_file SNP_MOTIF_FILE
                        SNP motif pair-wise list File, usually from step2
  -snp SNP_ID, --snp_id SNP_ID
                        SNP id

```
#### step3.3 SNP-TF radar plot  
**Example:** 
   ```
    # draw radar plot for instereing motif and SNP events
   MotifRaptor snpmotifradar \
   -wd step3_out \
   -sm step2_out/result_new_df_motifs_ENCFF512IML.txt \
   -sf step2_out/SNP_ENCFF512IML_features.txt \
   -pid rs7528684:MA0105.1__NFKB1
   ```

**Input:** Just specify the output directory from step2.
   
**Output:** Radar plot for the features/scores for interesting SNP-motif events.


<img src="https://github.com/pinellolab/MotifRaptor/blob/master/Document/pic5.png" alt="drawing" width="400"/>

**Usage:**
```
usage: MotifRaptor snpmotifradar [-h] [-wd WORKDIR] [-sm SNP_MOTIF_FILE]
                                 [-sf SNP_FEATURE_FILE] [-pid SNP_MOTIF_ID]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -sm SNP_MOTIF_FILE, --snp_motif_file SNP_MOTIF_FILE
                        SNP motif pair-wise list File, usually from step2
  -sf SNP_FEATURE_FILE, --snp_feature_file SNP_FEATURE_FILE
                        SNP feature file, usually from step2
  -pid SNP_MOTIF_ID, --snp_motif_id SNP_MOTIF_ID
                        SNP motif pair-wise ID
```


    

## Optional Step: Generating Database from a TF list.
Calculating disruption scores genome wide for a list of SNPs from a VCF file

<!--
### Install Scanner
1. Download files or git clone from this repository and keep them in the same folder
   ```
   git clone https://github.com/pinellolab/MotifRaptor.git
   
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
   *Optional*
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

 -->
