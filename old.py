__author__ = 'runsheng'
__data__ ="2015/07/14"
# The file from bam2coverage.py

# from bam file get the coverage as
# the following is the test code
import pysam
samfile = pysam.AlignmentFile("./4st_split_segment/cb12i_s.bam", "rb" )
for pileupcolumn in samfile.pileup("II",6029640,6029655):
    if pileupcolumn.pos in range(6029640,6029655+1):
        print ("\ncoverage at base %s = %s" %
                (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            print pileupread.query_position,
            print pileupread.indel
            print ('\tbase in read %s = %s' %
                    (pileupread.alignment.query_name,
                     pileupread.alignment.query_sequence[pileupread.query_position]))
samfile.close()


# The file bam2table.py
import pysam

# from bam file get the blast m6 like table
samfile = pysam.AlignmentFile("./4st_split_segment/cb4_bp_s.bam", "rb")

for read in samfile.fetch(): #the fetch is like to choose a subset of the bam file within a small region
    print read.query_name,'-' if read.is_reverse else '+',read.query_alignment_start,read.query_alignment_end,read.query_length,\
          samfile.getrname(read.reference_id), read.reference_start,read.reference_end,read.reference_length

samfile.close()
