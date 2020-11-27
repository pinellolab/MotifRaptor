# MotifRaptor

## Overview

##### Motivation: 
Genome-wide association studies (GWAS) have identified thousands of common trait-associated genetic variants but interpretation of their function remains challenging. These genetic variants can overlap the binding sites of transcription factors (TFs) and therefore could alter gene expression. However, we currently lack a systematic understanding on how this mechanism contributes to phenotype. 
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
## Configure Essential Databases

   Database contains essential data for general analysis, including DHS tracks, TF RNA-seq expressions, TF motifs, and TF pre-calucated scores.
   
   Please download the Database.zip from the links shown below. You can choose to download either a *testing database* or a *complete database*. The testing database contains all necessary files except that it only includes a handful number of TFs, in order for you to test the tutorial example on a regular machine. The complete database contains a complete list of TFs used in our real-world study, however, you need to have at least 220 G disk space, and we recommend running the programs on a cluster. 
   
   You may also choose to download the *testing database*, and then refer to the section *Build complete motif database* to compute a full TF list on your own, rather than downloading a big file. In this case, we also recommend running the programs on a cluster.
   
   
   ```
   # A testing database can be downloaded here.
   wget https://www.dropbox.com/s/gxeyzgl5m0u55w8/Database.zip
   unzip Database.zip
   ```
   
   ```
   # A complete database can be downloaded here. Make sure you have 220 G disk space.
   wget https://www.dropbox.com/s/kp5r82x55tfgawf/Database.zip
   unzip Database.zip
   ```
   
   Configuration for Motif-Raptor is to set up the absolute paths for general database and motif database.
   ```
   MotifRaptor set -pn databasedir -pv $PWD/Database/hg19/
   MotifRaptor set -pn motifdatabasedir -pv $PWD/Database/hg19/motifdatabase/
   ```
   
   Double check the paths are correctly set up.
   ```
   MotifRaptor info
   ```
   Skip the next section if you are only doing a tutorial.
   
## Build a Complete Motif Database
   
   (1) Download pfm files and change the motif database directory immediately.
   ```
   wget https://www.dropbox.com/s/hx9av7o16efxmus/motifdatabase.zip
   unzip motifdatabase.zip
   MotifRaptor set -pn motifdatabasedir -pv $PWD/motifdatabase
   ```
   *We use the JASPAR PFM format since it is well adopted and easy to download in batch from JASPAR database. This format is based on a simple flat file therefore it is easy to convert from/to. There are also computational tools that can be used directly to convert from other databases or file formats. For example, the biopython library has several parsers (Bio.motifs.pfm module) for PFM formats that can be used to generate a customized motif database.
   https://biopython.org/docs/1.75/api/Bio.motifs.pfm.html*
   
   (2) Download SNP list from 1000 Genome project
   ```
   wget https://www.dropbox.com/s/9gztf4mdblc44jo/1000G.EUR.QC.plink.simple.vcf
   ```
   The SNP list VCF file needs to have the first five columns ('CHR','POS','ID','REF','ALT') as follows:
   
   CHR | POS | ID | REF | ALT
   ------------ | ------------- | ------------- | ------------- | ------------- 
   1 | 2483961 | rs2258734 | A | G
   
   
   (3) Index the SNP list
   ```
   MotifRaptor snpindex -vcf 1000G.EUR.QC.plink.simple.vcf -gi genome_index -p 4
   ```
   
   (4) Scan motifs usnig pfm files on this SNP list
   ```
   MotifRaptor snpscan -gi genome_index -pfm ./motifdatabase/pfmfiles -mo ./motifdatabase/motifscanfiles -p 4
   ```
   *In the output folder, only '.scale' and '.score' files are useful. You may delete intermediate results in those folders.*
   ```
   cd motifdatabase/motifscanfiles/
   find ./ -type d -exec rm -rf '{}' \;
   ```


## Tutorial


### step 0. prepare input data (pre-processing) from GWAS summary statistics

**Input:** GWAS Summary Statistics 
Run Motif-Raptor from GWAS summary statistics. You may get summary statistics from UKBiobank, published paper, or other resources. These files may provide diffrent information. Please make sure the file contains the following columns. For the score column, Motif-Raptor currently supports pvalue, zscore, or chisquare.


ID | CHR | POS | REF | ALT | SCORE(pvalue, zscore, or chisquare)
------------ | ------------- | ------------- | ------------- | -------------  | ------------- 
rs2258734 | 1 | 2483961 | A | G | 0.003

**Example:**  Download the original data file from (Okada et al. 2010 Nature), and applying your own p-value cut-offs to define hits and nonhits. By default, p-value cutoff is 5E-8.
   This data file is ~450M. If your internet is limited, please download the zip file (~100M) and unzip it.
   ```
   wget https://www.dropbox.com/s/jnmpu63vqnlc0ig/RA_GWASmeta_TransEthnic_v2.txt
   
   #alternative zip file
   wget https://www.dropbox.com/s/c194x1z0bhntfbs/RA_GWASmeta_TransEthnic_v2.zip
   unzip RA_GWASmeta_TransEthnic_v2.zip
   ```
   <!--

   ```
   #original link of this data file (you don't need to download it again):
   #wget https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Okada/RA_GWASmeta_TransEthnic_v2.txt.gz
   #gunzip RA_GWASmeta_TransEthnic_v2.txt.gz
   ```

   -->
   
   In this file, columns 1,2,3,4,5,9 are ID,CHR,POS,REF,ALT,SCORE as defined above. Here the score is pvalue.
   ```
   MotifRaptor preprocess -gwas RA_GWASmeta_TransEthnic_v2.txt -cn 1,2,3,4,5,9 -st pvalue -th 5E-8
   ```
   

**Output:** Information for SNP hits and non-hits. *hitSNP_list.txt*  *nonhitSNP_list.txt* and *hitSNP_list.vcf*
For the example from (Okada et al. 2010 Nature), you can also download our processed results, if you haven't run on your own.
   ```
   wget https://www.dropbox.com/s/gpnudp1ba4d2gq3/hitSNP_list.txt
   wget https://www.dropbox.com/s/7dfrph1dnrad894/nonhitSNP_list.txt
   wget https://www.dropbox.com/s/73seqi42hgodupg/hitSNP_list.vcf
   ```
**Double Check Output Format:** *hitSNP_list.txt* and *nonhitSNP_list.txt* are two files with the following format:

ID | CHR | POS
------------ | ------------- | ------------- 
rs2258734 | 1 | 2483961

*hitSNP_list.vcf* is the file with the following format (with two more columns of the polymorphism information:

ID | CHR | POS | REF | ALT
------------ | ------------- | ------------- | ------------- | ------------- 
rs2258734 | 1 | 2483961 | A | G

**Usage:**
   ```
   usage: MotifRaptor preprocess [-h] [-gwas SUMMARYSTATSFILE]
                              [-cn COLUMN_NUMBERS] [-st SCORE_TYPE]
                              [-th SCORE_PVALUE_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -gwas SUMMARYSTATSFILE, --gwas_summary SUMMARYSTATSFILE
                        GWAS summary statistic file
  -cn COLUMN_NUMBERS, --column_numbers COLUMN_NUMBERS
                        provide six column numbers for information in such
                        order: ID,CHR,POS,A1,A2,SCORE eg. 1,2,3,4,5,6
  -st SCORE_TYPE, --score_type SCORE_TYPE
                        Score type in GWAS summary statistic file: pvalue or
                        zscore or chisquare
  -th SCORE_PVALUE_THRESHOLD, --threshold SCORE_PVALUE_THRESHOLD
                        threads for pvalue - default 5E-8
   ```

**Bonus function:** 

Motif-Raptor allows user to preprocess from the UKBB (v3) summary statistics TSV files directly using a much simpler command.

   ```
   wget https://www.dropbox.com/s/axthfv12j7pbav4/30690_raw.gwas.imputed_v3.both_sexes.tsv
   MotifRaptor preprocess_ukbb_v3 -gwas 30690_raw.gwas.imputed_v3.both_sexes.tsv -th 5E-8
   ```
   
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

## Conclusion

In summary, Motif-Raptor is a computational toolkit to test the significance of the effects of genetic variants from GWAS analyses on transcription factor binding sites. We believe that its adoption will help the genomic community in prioritizing potential cell type-specific, causal variants from GWAS summary statistics and to generate important hypotheses and insights to the mechanisms of action of genetic variants in complex disease.

## Contact

Luca Pinello: lpinello at mgh.harvard.edu

Qiuming Yao: yao.ornl at gmail.com

