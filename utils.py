#filename: utils.py
#prepare some functions that can be used globely, such as read the reference to the memory, read the fastq file to dict

#get the fastq file into dict
from Bio import SeqIO
def fastq2dic(fastqfile):
    """
    Give a fastq file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large readfile may take a lot of RAM
    """
    handle = open(fastqfile, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
    handle.close()
    return record_dict

def fasta2dic(fastafile):
    """
    Give a fasta file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large genome may take a lot of RAM
    """
    handle = open(fastafile, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
    handle.close()
    return record_dict

def chr_select(record_dict, chr, start,end):
    """
    Note the start and end is 0 based
    give the name of refdic, and the chr, start and end to be used
    return the name and sequence
    for example, chrcut(record_dict, "I", 0,100) returns
    """
    name=chr+ ":"+str(start)+"_"+str(end)
    seq=str(record_dict[chr].seq)[start:end]
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

import subprocess, signal, logging, os

class Alarm(Exception):
    pass
def alarm_handler(signum, frame):
    raise Alarm

def exe(cmd, timeout=-1):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtools orders
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, close_fds=True,
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout))
    try:
        stdoutVal, stderrVal =  proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d" % (timeout)))
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return 214,None,None
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

def bwaindex(fastafile):
    pass

def bwamem():
    cmd="bwa mem"

### unit test code ###
if __name__=="__main__":
    print "aa"