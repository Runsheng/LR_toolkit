# filename: check_N.py
# restore the database
# post-assembly function
"""
From the pickle file ,restore the N_list

for a position of N_list, get the reads overlap the region, with =-10 flanking region

and check the Ns within a length range,
give the position and the reads that can possibly overlapped the region
take the reads from the fatsq file and prepare them to a simple assembler, so the gap can be filled
"""
import cPickle as pickle
import pysam
from utils import fasta2dic

# unused function
def get_readname(samfile_name, chr, start,end):
    """
    :param samfile_name:the name of samfile to be parsed
    :param chr: str
    :param start: int
    :param end: int
    :return:a list contain the names of the reads overlap the region
    """
    read_names=[]
    with pysam.AlignmentFile(samfile_name) as samfile:
        for read in samfile.fetch(chr,start-10,end+10):
            read_names.append(read.query_name)
    return read_names

# unused function
def write_fasta(read_names, read_dic, out="support_read.fasta"):
    """
    The function to write selected reads into fasta file
    input: a list containing all the read names to be written,
           a dictionary containing all the {readname:sequence} in biopython SeqIO format
    output:a file
    Note: the biopython SeqIO is really toooo complex to be used in writing an simple fasta file
    """
    with open(out,"w") as f:
        for name in read_names:
            sequence = str(read_dic[name].seq)
            f.write(">")
            f.write(name)
            f.write("\n")
            f.write(sequence)
            f.write("\n")

# post-run function
def main():
    with open("./4st_split_segment/N_list.dat","rb") as fp:
        N_list = pickle.load(fp)
        for i in range(1000,2000):
            if 10000<N_list[i][2]-N_list[i][1]:
                print N_list[i]
                chr, start, end = N_list[i]
                with pysam.AlignmentFile("./4st_split_segment/cb12i_s.bam", "rb") as samfile:
                    for read in samfile.fetch(chr,start-10,end+10):
                        print read.query_name
    print "The N_list contains", len(N_list), "gaps."


if __name__ == "__main__":
    main()