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

##### At this step, only pysam is used
- 1. Summary_N, get all the Nsite, store as a list, pickle dump as a file, N.dat, 
- 2. cons, using the mapped bam, scan the N region, get new cons of these region, 
- 3. fillN, fill the N with the cons, 
- 4. Summary_N of newly filled genome, to get a N_list again, 
        compare the N_list with the former one, 
        - if identical (or almost identical, for example, change <1000 bp), stop
        - if not, repeat 1-4
