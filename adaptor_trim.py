# internal package import
from utils import fasta2dic, fastq2dic, timer
from ssw.ssw_wrap import Aligner

# third part import
import pysam
from Bio import SeqIO

# use the adaptor library inside code? -> hardcore first, then sth else
>TruSeq_Universal_Adapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_Index_Adapter
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC


def _write_fastq(fastq_outpath,start,end):
    w_hdl = open(fastq_outpath, 'w')
    w_seq = seq[start:end]
    SeqIO.write([w_seq], w_hdl, 'fastq')
    w_hdl.close()


def ssw_ref(ref_path, read_path):
    ref_dict=fasta2dic(ref_path)
    read_dict=fastq2dic(ref_path)
    for ref_name, ref_read in ref_dict.iteritems():
        for

seq1="TTGAATTCCTGCCAAAAAAATTTTCTTCAAAAATTTTAATTTCCAGCCAAAATGTTTTTTTTCCGAAAATTTAAGTTTCCCGCTAAAATATTTTTCTCAAAAATTTTAATTTTCCGCCAAAAAATTTTCCAGAAAGTTTGAATTTCTTCAAAAAAT"
seq2="ATTTTCCAGAAAGTTTGAATTTCTTCAAAAAATTTGAATTCCTGCCAAAAAAATTTTCTTCAAAAATTTTAATTTCCAGCCAAAATTTTTTTTTCGAAAATTTAAATTTCCCGCTAAAATATTTTTCTCAAAAATTTTAATTTCCCGCCAAAAATATT"
aligner = Aligner(seq1, report_cigar=True)
aln=aligner.align(seq2)

print aln.ref_begin,aln.ref_end,aln.query_begin,aln.query_end,aln.cigar_string


# make the alignment first

if __name__=="__main":
    pass
    print("at least this runs")

