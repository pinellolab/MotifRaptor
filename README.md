# MotifRaptor

## Demo
This is an example of characterizing interesting TFs for Rheumatoid Arthritus GWAS data.
Follow the steps in Demo/RA/RA_demo.ipynb

## Installation
1. Download files from this repository and keep them in the same folder

2. Extra python packages:
   ```
    pip install twobitreader
   ```
3. Install mksary

   3.1 Download [libdivsufsort](https://goo.gl/hUjvMF), which is the library available [in GitHub](https://github.com/y-256/libdivsufsort) changed in the file "mksary.c" which now includes the computation of the LCP array.
   
   3.2 Build   
   ```
    cd libdivsufsort 
    mkdir build 
    cd build
    make 
    sudo make install
    ```

    3.3 Copy the executable "mksary" into your python-working directory

4. Compile the cython code: cythonize -a -i motif_matching_lcp.pyx
## Data
1. VCF file:

   Download VCF file from [here](https://www.dropbox.com/s/9gztf4mdblc44jo/1000G.EUR.QC.plink.simple.vcf?dl=0)

2. Genome 2Bit File:
   
   Download from UCSC website [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) 

## Building disruption scores genome wide for a list of SNPs from a VCF file

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
