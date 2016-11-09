#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 11/7/16 2:39 PM
# @Author  : Runsheng
# @File    : muti_mpileup_test.py

from utils import myexe, dic2fasta
import vcf

from logger import log_summary
logger=log_summary()


def split_genome(ref_dict, interval=500000):
    """
    Get a normal intervals with length=interval, to feed into mpileup
    Note: VCF is 0 based file, but SAM is 1 base file
    :param ref_dict:
    :param interval:
    :return: as chr:start-end
    """
    genome_pieces=[]

    for k, v in ref_dict.items():
        if len(v)>interval:
            for n in range(len(v)//interval+1):
                if interval*(n+1)<=len(v):
                    p_str="{chro}:{start}-{end}".format(chro=k,start=interval*n+1,end=interval*(n+1)+1)
                else:
                    p_str="{chro}:{start}-{end}".format(chro=k,start=interval*n,end=len(v))
                genome_pieces.append(p_str)
        else:
            p_str="{chro}:{start}-{end}".format(chro=k,start=1,end=len(v))
            genome_pieces.append(p_str)

    return genome_pieces


def _test_caller_single(wkdir, ref,bamfile,region):
    pass


def _parser_upper(strings):
    """
    used only for the insertions call using a lower cased genome
    the insertion will showed as UPPER in the record.ALT
    """
    upper_l=[]
    for i in strings:
        if i.isupper():
            upper_l.append(i)
    return "".join(upper_l)


def get_ins_region(vcf_file, len_cutoff=5):
    """
    The filter is built inside, with
    1. homo deletion
    2. coverage? DP>?
    Normal vcf is 0 based
    :return: a dict as {I:100^101: "TTAC"}
    """
    bed_d={}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.var_subtype=="ins" and record.aaf[0]==1:  # choose only the homo insetion
            alt_seq=str(record.ALT[0]).upper()
            ref_seq=str(record.REF).upper()
            insertion=alt_seq[len(ref_seq):]

            #print alt_seq, ref_seq, insertion
            ins_len=len(insertion)
            if ins_len>=len_cutoff:
                ins_start=record.affected_start
                name=record.CHROM+":"+str(ins_start)+"^"+str(ins_start+1)
                bed_d[name]=insertion
    return bed_d


def get_ins_region_clc(vcf_file, len_cutoff=5):
    """
    The filter is built inside, with
    1. homo deletion
    2. coverage? DP>?
    BCFtools vcf is 1 based, CLC is 0 based?
    :return: a dict as {I:100^101: "TTAC"}
    """
    bed_d={}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.var_subtype=="ins" and record.aaf[0]==1:  # choose only the homo insetion
            alt_seq=str(record.ALT[0]).upper()
            ref_seq=str(record.REF).upper()
            insertion=alt_seq[len(ref_seq):]

            #print alt_seq, ref_seq, insertion
            ins_len=len(insertion)
            if ins_len>=len_cutoff:
                ins_start=record.affected_start
                name=record.CHROM+":"+str(ins_start+1)+"^"+str(ins_start+2)
                bed_d[name]=insertion
    return bed_d


def _read_insertion(filename):
    """
    Only used to read a bed-like insertion file and store it with position in key
    used in debug using the other pipeline out of

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


def write_insertion_used(ins_d, outfile):
    with open(outfile, "w") as fw:
        for k,v in ins_d.items():
            chro=k.split(":")[0]
            start = int(k.split(":")[1].split("^")[0])
            end= start+1
            line=[chro, str(start),str(end),v]
            fw.write("\t".join(line))
            fw.write("\n")


def write_insertion_fasta(record_dict, insertion, outfile="inserted.fasta"):
    new_dict={}
    for name in record_dict.keys():
        subinsertion = {}  # note , a dict is not really needed, a list can do the job
        for key in insertion.keys():
            if name == key.split(":")[0]:
                # use number as index, perform better in sorting
                start_key = int(key.split(":")[1].split("^")[0])
                subinsertion[start_key] = insertion[key]
        keys_in_order = sorted(subinsertion.keys())
        seq = list(str(record_dict[name].seq))
        i = 0  # i count for how many insertions were added
        for key in keys_in_order:
            seq.insert((key - 1 + i), subinsertion[key])
            i += 1
        if i!=0:
            print("%d insertions are filled to %s" %(i, name))

        new_dict[name]="".join(seq)
    dic2fasta(new_dict, outfile)

if __name__=="__main__":
    ref="/home/zhaolab1/myapp/LR_toolkit/test/wkdir/temp/ref/round7_nfill.fasta"
    bamfile="/home/zhaolab1/myapp/LR_toolkit/test/wkdir/temp/round8_s.bam"
    wkdir="/home/zhaolab1/myapp/LR_toolkit/test/wkdir/"

    from utils import fasta2dic, my_chdir

    import os
    os.chdir(wkdir)

    ref_dict=fasta2dic("/home/zhaolab1/myapp/LR_toolkit/test/wkdir/temp/ref/round7_nfill.fasta")
    genome_pieces=split_genome(ref_dict, 500000)
    #from pprint import pprint
    #pprint(genome_pieces)

    #######set up
    vcf_dir = os.path.join(wkdir, "temp/vcf")
    #os.makedirs(vcf_dir)
    os.chdir(vcf_dir)


    ###### run code
    #def wrapper_bam2vcf_map(region):
    #    wrapper_bam2vcf_single(ref, bamfile, region, prefix="default")

    #pool=Pool(32)
    #pool.map(wrapper_bam2vcf_map, genome_pieces)


    ##### filter result and clean up
    from glob import glob
    vcf_files=glob("*.vcf")
    import vcf
    ins_all={}
    for vcf_one in vcf_files:
        ins_one=get_ins_region(vcf_one, len_cutoff=5)
        ins_all.update(ins_one)
    write_insertion_used(ins_all, "ins_round8.txt")
    write_insertion_fasta(ref_dict,ins_all,"round8_inserted.fasta")