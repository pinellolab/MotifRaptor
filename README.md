# MotifRaptor
1. Compile the cython code: cythonize -a -i motif_matching_lcp.pyx
2. Run the help to test: python Motif_Scan.py -h
3. Two modes of command in Motif_Scan: index, pfmscan
3.1 Example 1. Index a genome: python Motif_Scan.py index -g genome_database -p number_of_threads
3.2 Example 2. Scan motif PFM files: python Motif_Scan.py pfmscan -gi indexed_genome_database -pfm motif_pfm_folder -mo motifscan_result -p number_of_threads

