# filename: utils.py
# prepare some functions that can be used globally, such as read the reference to the memory, read the fastq file to dict
import subprocess
import signal  # only used for the exe function
import time
from functools import wraps

import os

# The sequence operation functions--------------------------------------------------------------------------------------
from Bio import SeqIO
from Bio import SeqIO
import gzip

def fastq2dic(fastqfile):
    """
    Give a fastq file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large readfile may take a lot of RAM
    """
    if ".gz" in fastqfile:
        handle=gzip.open(fastqfile, "rU")
    else:
        handle=open(fastqfile, "rU")
    record_dict=SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
    handle.close()
    return record_dict


def fasta2dic(fastafile):
    """
    Give a fasta file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large genome may take a lot of RAM
    """
    handle=open(fastafile, "rU")
    record_dict=SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
    handle.close()
    return record_dict


def dic2dic(record_dict):
    """
    :param record_dict: a SeqIO dict generated by Biopython
    :return the dict contain {name:seq}
    """
    seq_dict={}
    for k,v in record_dict.iteritems():
        seq=str(v.seq)
        seq_dict[k]=seq
    return seq_dict


def chr_select(seq_dict, chro, start,end):
    """
    Note the start and end is 0 based
    give the name of refdic, and the chr, start and end to be used
    return the name and sequence
    for example, chrcut(record_dict, "I", 100,109) returns
     ("I:100_109","AAAAAAAAAA")
    """
    name=chro+ ":"+str(start)+"_"+str(end)
    seq=str(seq_dict[chro][start:end].seq)
    return name,seq


def dic2fasta(record_dict,out="record_dict.fasta"):
    with open(out,"w") as f:
        for record in record_dict.keys():
            name=record
            seq=record_dict[name]
            f.write(">")
            f.write(name)
            f.write("\n")
            f.write(seq)
            f.write("\n")

def max_get(listname):
    """
    :param listname, a list containning some sequence strings;
    :return: the longest sequence string in the list.
    """
    chose=None
    max_length=len(listname[0])
    for n in range(1,len(listname)):
        if len(listname[n])>max_length:
            max_length=len(listname[n])
    for element in listname:
        if len(element)==max_length:
            chose=element
    return chose
# ----------------------------------------------------------------------------------------------------------------------

# The system operation functions----------------------------------------------------------------------------------------

def myexe(cmd, timeout=0):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtool orders
    """
    def setupAlarm():
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(timeout)

    def alarmHandler(signum, frame):
        sys.exit(1)

    proc=subprocess.Popen(cmd, shell=True, preexec_fn=setupAlarm,
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err=proc.communicate()
    print err
    return out, err, proc.returncode


def timer(fn):
    """
    Used for debug
    :param fn: the function
    :print: the time the function used
    """
    @wraps(fn)
    def wrapper(*args, **kwargs):
        ts = time.time()
        result = fn(*args, **kwargs)
        te = time.time()
        print "function = {0}".format(fn.__name__)
        print "    time = %.6f sec" % (te-ts)
        return result
    return wrapper
# ----------------------------------------------------------------------------------------------------------------------


def bwa_index(ref):
    myexe("bwa index %s" % ref)

def bwa_mem(ref, reads, core=15):
    '''
    :param core: the core used to run bwa mem
    :param ref: the full path or relative path to the ref file, need to be in fasta format
    :param reads: the full path or relative path to the read file, need to be in fastq format
    :return: None? or sth?
    '''
    assert ref.split(".")[-1] in ["fa","fasta"], "Check your reference type, and make them end with .fa or .fasta"
    assert reads.split(".")[-1] in ["fq","fastq"], "Check your reads type, and make them end with .fq or .fastq"

    ref_short=ref.split("/")[-1].split(".")[0]
    reads_short=reads.split("/")[-1].split(".")[0]
    reads_short_s=reads_short+"s"

    # index
    myexe("bwa index %s" % ref)
    # bwa to bam file
    print myexe("bwa mem -t {core} {refpath} {readpath} | samtools view -bS - |samtools sort - {reads_short_s}."
          .format(core=core, readpath=reads, refpath=ref, reads_short_s=reads_short_s))
    # sort bam file
    # myexe("samtools sort {reads_short}.bam {reads_short_s}".
    #      format(reads_short=reads_short,reads_short_s=reads_short_s))
    # index bam file
    myexe("samtools index {reads_short_s}.bam".
          format(reads_short_s=reads_short_s))
    # clean up
    #print myexe("rm {reads_short}.bam".format(reads_short=reads_short))



### unit test code ###
if __name__== "__main__":
    print myexe("ls")[0]