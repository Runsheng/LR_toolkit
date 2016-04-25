
from utils import myexe
import gzip
from Bio import SeqIO

# the 'SRR003265.filt.fastq.gz' shoule be in the path, or else get it with
# myexe("wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz")

# for the ones without zip, just read and subset them

def main():
    recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz'), 'fastq')
    fw=open("aa.fastq","w")
    for rec in recs:
        # select 10 nucleotides, and write them into a single file/ one by one
        SeqIO.write(rec[4:14], fw, 'fastq')
    fw.close()

if __name__=="__main__":
    # main()
    # use the following sentence to write a adaptor fasta file
    with open("adaptor.fasta","w") as f:
        f.write(">TruSeq_Universal_Adapter\n")
        f.write("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n")
        f.write(">TruSeq_Index_Adapter\n")
        f.write("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n")

