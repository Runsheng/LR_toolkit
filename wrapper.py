#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/28/16 3:33 PM
# @Author  : Runsheng     
# @File    : wrapper.py

from utils import myexe
from logger import log_summary
from functools import wraps

logger=log_summary()


def show_runcmd(fn):
    """
    used for showing the running commander in myexe
    :param fn:
    :return:
    """
    @wraps(fn)
    def wrapper(*args,**kwargs):
        logger.info("RUNNING    "+str(" ".join(*args)))
        results=fn(*args,**kwargs)
        return results
    return wrapper



def wrapper_bwamem(index, read_list, prefix="default", core=32, w=15000, k=50):
    """
    general wrapper for bwa mem to be used in python cmd line
    para: index, the bwa index of reference fasta file
    para: read_list, the fastq file names in a python list
    para: prefix, the prefix for the output bam file
    para: core, the cores used for mapping and sorting

    return: None
    """
    prefix=read_list[0].split("_")[0] if prefix=="default" else prefix
    read_str=" ".join(read_list)
    print myexe("pwd")
    myexe("bwa mem -t {core} -w {w} -k {k} {index} {read_str} > {prefix}.sam".format(core=core, w=w, k=k,
                                                                                      index=index, read_str=read_str,prefix=prefix))
    myexe("samtools view -bS {prefix}.sam >{prefix}.bam".format(prefix=prefix))
    myexe("samtools sort -@ {core} {prefix}.bam -o {prefix}_s.bam".format(core=core,prefix=prefix)) # only valid for samtools >1.2
    myexe("samtools index {prefix}_s.bam".format(prefix=prefix))
    myexe("rm {prefix}.sam {prefix}.bam".format(prefix=prefix))


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



def __wrapper_bam2vcf(ref, bamfile, prefix="default", core=32):
    """
    para: ref, the reference fasta file
    para: bamfile: a single sorted bam file, with index
    para: prefix, the prefix for the output vcf file

    return: None
    """
    print myexe("pwd")
    prefix=bamfile.split(".")[0].split("_")[0] if prefix=="default" else prefix
    cmd_mpileup=("samtools mpileup -ugf {ref} {bam} | bcftools call -vmO z -o {prefix}.bcf"
                 .format(ref=ref,bam=bamfile,prefix=prefix))
    cmd_tovcf=("bcftools view {prefix}.bcf>{prefix}.vcf".
               format(prefix=prefix))

    logger.info("RUNNING  "+cmd_mpileup)
    myexe(cmd_mpileup)
    logger.info("RUNNING  "+cmd_tovcf)
    myexe(cmd_tovcf)
    print "++++++++done with VCF calling++++++++++"



def __bwa_index(ref):
    myexe("bwa index %s" % ref)

def __bwa_mem(ref, reads, core=15):
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
