# To add the insertion sequence into the reference sequence
# the bed file contains the "chr\tstart\tend\tsequence\n" chrs
# the fasta file is the reference sequence to be filled
# require Biopython:SeqIO to parse the fasta file
# runsheng 2014/12/15
# insert new strings into old strings, make a list to store the number of the sequence 
from utils import fasta2dic

def read_insertion(filename):
    """
    :param filename is the bed file name for the insertion, as "chr\tstart\tend\tsequence\n"
    :return a dict as {insertion:seq}, for example {"I:130^131":"AAATTTTCCC"}
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

# usage todo: need to be modified to be main function
insertion=read_insertion("insertion_filled.txt")
record_dict=fasta2dic("cb4_nfilled.fasta")

def write_insertion(outfile="inserted.fasta"):
    f=open(outfile)
    for name in record_dict.keys():
        subinsertion={}
        print name
        for key in insertion.keys():
            if name==key.split(":")[0]:
                subinsertion[key]=insertion[key]
        seq_old=str(record_dict[name].seq)
        seq_new=[]
        for i in range(0,len(seq_old)):
            for key in subinsertion.keys():
                if i==(int(key.split("^")[1])-2):
                    seq_new.append(subinsertion[key])
            seq_new.append(seq_old[i])
        print "old", len(seq_old)
        print "new", len(seq_new)
        f.write(">"+str(name)+"\n")
        f.write("".join(seq_new))
        f.write("\n")
    f.close() b
