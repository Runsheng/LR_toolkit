# internal package import
from utils import fasta2dic, timer, bwa_mem,myexe
from ssw.ssw_wrap import Aligner

# third part import
import pysam
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from cStringIO import StringIO
from Bio.Blast import NCBIXML

# the adaptor reads have been written into adaptor.fasta file
# for the illumina truseq reads, only the universal(5') and adaptor(3') is enough,
# and the index for the 384 barcoding is single-ended, with 8 nucleotides

def _ssw_ref(read_path, ref_path="adaptor.fasta"):
    '''
    The orign method to find the adaptor position, slow and buggy
    :param read_path:
    :param ref_path:
    :return:
    '''
    ref_dict=fasta2dic(ref_path)
    fqs_in=SeqIO.parse(open(read_path),'fastq')
    fqs_out=open(read_path.replace(".fq","").replace(".fastq","")+"trimmed"+".fastq","w")
    i=0
    for fq in fqs_in:
        for ref_name, ref_read in ref_dict.iteritems():
            len_cutoff=len(ref_read)-3
            aligner=Aligner(ref_read,report_cigar=True)
            aln=aligner.align(fq, min_score=len_cutoff*2-10, min_len=len_cutoff)
            if aln !=None:
                print ref_name, aln.score, aln.ref_begin,aln.ref_end,aln.query_begin,aln.query_end,aln.cigar_string
    fqs_out.close()
    # run with min_score=len_cutoff*2-10, finished in 80 min but found no match


def bamview(bamfile):
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for read in samfile.fetch():
            print read.cigarstring


def adaptor_blast(query,dbpatch="adaptor.fasta"):
    # build the blast db, maybe adding an asserting to identify the exsentise of the db is better
    db=dbpatch.split(".")[0]
    print myexe("makeblastdb -in %s -dbtype nucl -input_type fasta -out %s" % (dbpatch,db))

    blastn_cline = NcbiblastnCommandline(db=db, outfmt=5)
    out, err = blastn_cline(stdin=query)
    blast_records = NCBIXML.read(StringIO(out))  # return is a generator, need a loop to parse the result
    return blast_records



# make the alignment first

if __name__ == "__main__":
    print("Running trim process")
    # ssw test code
    # ssw_ref("cb12.fq","adaptor.fasta")

    # bwa test code
    # bwa_mem(ref="adaptor.fasta", reads="cb12.fq", core=15)
    # print bamview("cb12_s.bam")

    # blast test code
    # for the determination
    # for the detection of the contig sinside the