#from bam file get the blast m6 like table
import pysam
samfile = pysam.AlignmentFile("cb4_bp_s.bam", "rb")

for read in samfile.fetch(): #the fetch is like to choose a subset of the bam file within a small region
    print read.query_name,'-' if read.is_reverse else '+',read.query_alignment_start,read.query_alignment_end,read.query_length,\
          samfile.getrname(read.reference_id), read.reference_start,read.reference_end,read.reference_length

samfile.close()


# to get per base coverage at a place
# the str of the pileup function
# Class: AlignmentFile.pileup.pileups
# Classname: AlignmentFile.pileup==pileupcolumn; AlignmentFile.pileup.pileups==pileupread
# Note: the N_list do not include the end of the gap region, so is needed to be expanded

#todo: change the marix in this one
import pysam
print N_list[190]
chro, start,end =N_list[190]
samfile = pysam.AlignmentFile("cb12i_s.bam", "rb")
print chr_select(record_dict, chro,start,end+1)

for pileupcolumn in samfile.pileup(chro,start,end+1, truncate=True):
    print ("\ncoverage at base %s = %s" %
            (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        print pileupread.alignment.query_sequence[pileupread.query_position]
        print pileupread.is_del  #if this position is a deletion,then the replacemnet can be ignored
        print pileupread.indel
samfile.close()