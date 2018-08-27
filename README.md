# MotifRaptor

## Installation
1. pip install twobitreader
2. Install mksary

   2.1 Download [libdivsufsort](https://goo.gl/hUjvMF)
   
   2.2 Build   
   ```
    cd libdivsufsort 
    mkdir build 
    cd build
    make 
    sudo make install
    ```

    2.3 Copy the executable "mksary" into your python-working directory

3. Compile the cython code: cythonize -a -i motif_matching_lcp.pyx
## Data
1. VCF file:

   Download VCF file from [here](https://www.dropbox.com/s/9gztf4mdblc44jo/1000G.EUR.QC.plink.simple.vcf?dl=0)

2. Genome 2Bit File:
   
   Download from UCSC website [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) 

## Test run
1. Run the help to test: python Motif_Scan.py -h

2. Pre-process the VCF file:
     ```
      python Process_vcf.py
     ```

3. Two modes of command in Motif_Scan: index, pfmscan

    3.1 Example 1. Index a genome: 

     ```
      python Motif_Scan.py index -g genome_database -p number_of_threads
     ```
    
    3.2 Example 2. Scan motif PFM files: 

     ```
      python Motif_Scan.py pfmscan -gi indexed_genome_database -pfm motif_pfm_folder -mo motifscan_result -p number_of_threads
     ```
