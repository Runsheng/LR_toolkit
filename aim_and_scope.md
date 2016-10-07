##Aim for this LR_toolkit
- Input: accurate long reads, for example 
    - illumina synthetic long read or 
    - corrected pacbio reads
- Input: current genome assembly 
    - fill the gap
    - list some re-arrangements and correct them 

- Output: improved genome
    - fasta file
    - gff file and summary of the fixed region
        - The gff should include the newly added region
        - the intersection between the added sequence and raw gff file need to be generated
        - maybe a chain file could be useful

###The flow for the pipeline is:

- Summary_N, get all the Nsite, store as a file, pass to cons  
