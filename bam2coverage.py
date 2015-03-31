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