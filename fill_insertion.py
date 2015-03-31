# To add the insertion sequence into the reference sequence
# the bed file contains the "chr\tstart\tend\tsequence\n" chrs
# the fasta file is the reference sequence to be filled
# require Biopython:SeqIO to parse the fasta file
# runsheng 2014/12/15
# insert new strings into old strings, make a list to store the number of the sequence 
from utils import reftodic


def read_insertion(filename):
    """
    read the insertion table in this bed-like format
    chr\tstart\tend\tsequence\n
    into a dict
    """
    d={}
    file_open=open(filename)
    lines=file_open.readlines()
    for line in lines:
        # store the name in the format of key:value
        # {"I:130^131":"AAATTTTCCC"}
        insertion=line.split("\t")[0]+":"+line.split("\t")[1]+"^"+line.split("\t")[2]
        sequence=line.split("\t")[3].strip()
        d[insertion]=sequence
    file_open.close()
    return d

# usage
insertion=read_insertion("insertion_filled.txt")
record_dict=reftodic("cb4_nfilled.fasta")

def write_insertion(outfile="inserted.fasta"):
    f=open(outfile)
    for name in record_dict.keys():
        subinsertion={}  # note , a dict is not really needed, a list can do the job
        print name
        for key in insertion.keys():
            if name==key.split(":")[0]:
                start_key=int(key.split(":")[1].split("^")[0]) # use number as index, perform better in sorting
                subinsertion[start_key]=insertion[key]
        keys_in_order=sorted(subinsertion.keys())
        seq=list(str(record_dict[name].seq))
        i=0 # i count for how many insertions were added
        for key in keys_in_order:
            seq.insert((key-1+i),subinsertion[key])
            i+=1
        f.write(">"+str(name)+"\n")
        f.write("".join(seq))
        f.write("\n")
    f.close()
