import pysam

# from bam file get the blast m6 like table
samfile = pysam.AlignmentFile("./4st_split_segment/cb4_bp_s.bam", "rb")

for read in samfile.fetch(): #the fetch is like to choose a subset of the bam file within a small region
    print read.query_name,'-' if read.is_reverse else '+',read.query_alignment_start,read.query_alignment_end,read.query_length,\
          samfile.getrname(read.reference_id), read.reference_start,read.reference_end,read.reference_length

samfile.close()


