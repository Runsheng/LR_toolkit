
import pysam

# samfile=pysam.AlignmentFile("cb12nin_s.bam", "rb")

def edge_parse(samfile):
    """
    Input: the samfile is a pysam.Alignmentfile object
    Output: the bp_list liked tuple-list
    """
    ref_length=dict(zip(samfile.references,samfile.lengths))

    for ref in ref_length.keys():
        clipping_l=[]
        clipping_r=[]
        for read in samfile.fetch(ref, 0, 100):
            # cigartuples, 4 is softclipping, 5 is hardclipping,
            # TODO: the hard clipping has been ignored, can only be recovered when the fastq file directly
            if read.cigartuples[0][0]==4:
                clipping_l.append(read.seq[0:int(read.cigartuples[0][1])])
        for read in samfile.fetch(ref, (ref_length[ref]-100),ref_length[ref]):
            if read.cigartuples[-1][0]==4:
                clipping_r.append(read.seq[0:int(read.cigartuples[0][1])])
        try:
            max_get(clipping_l) # to avoid the output of empty sequence name line
        except Exception as e:
            print e
            print "No left clipping signal found in chr %s ." % ref
        try:
            max_get(clipping_r)
        except Exception as e:
            print e
            print "No right clipping signal found in chr %s ." % ref

if __name__=="__main__":
    samfile=pysam.AlignmentFile("cb12nin_s.bam", "rb")
