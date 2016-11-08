#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 11/7/16 2:39 PM
# @Author  : Runsheng
# @File    : muti_mpileup_test.py

from utils import myexe
from multiprocessing import Pool

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


def wrapper_bam2vcf_single(ref, bamfile, region, prefix="default"):
    """
    get a single vcf caller wrapper that can easy be called using multiprocess
    :param ref:
    :param bamfile:
    :param region:
    :param prefix:
    :return:
    """
    print(myexe("pwd"))
    prefix=bamfile.split("/")[-1].split("_")[0] if prefix=="default" else prefix
    print(prefix)
    outname=prefix+"_"+region
    cmd_mpileup=("samtools mpileup -r {region} -ugf {ref} {bam} | bcftools call -v -m -O v -o {prefix}.vcf"
                 .format(region=region, ref=ref,bam=bamfile,prefix=outname))

    logger.info("RUNNING  "+cmd_mpileup)
    myexe(cmd_mpileup)
    logger.info("++++++++done with VCF calling++++++++++")


def test_caller_single(wkdir, ref,bamfile,region):
    pass


def _parser_upper(strings):
    """
    used only for the insertions,
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
    :return: a dict as {I:100^101: "TTAC"}
    """
    bed_d={}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.var_subtype=="ins" and record.aaf[0]==1:  # choose only the homo insetion
            alt_seq=record.ALT[0]
            insertion=parser_upper(alt_seq)

            ins_len=len(insertion)
            if ins_len>=5:
                ins_start=record.affected_start
                name=record.CHROM+":"+str(ins_start)+"^"+str(ins_start+1)
                bed_d[name]=insertion
    return bed_d



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
    os.makedirs(vcf_dir)
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
    for vcf_one in vcf_files:
        vcf_reader=vcf.Reader(open(vcf_one), "r")


