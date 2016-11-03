#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/28/16 10:22 AM
# @Author  : Runsheng
# @File    : test_mapping.py


# define the main flow of the test run for cb4 genome
# firstly do the mapping
from __future__ import print_function
import os
import shutil
import logging
import glob

# used in run_mapper
from utils import myexe
from wrapper import wrapper_bwamem
from wrapper import wrapper_bam2vcf
from utils import myglob  # recursive glob for python 2
from glob import glob # normal

# used in run_nreplace
from summary_N import summary_N
from utils import fasta2dic
from cons import write_nreplace

from sequence_replace import read_nreplace,sequence_replace, write_nreplace_used


# debug and profiling functions
from utils import timer
from wrapper import show_runcmd

# ----------------------------------------------------------------------------------------------------------------------
def pre_dir_file(i,work_dir,ref_file=None):
    """
    prepare the dir and files for the i round mapping and fill
    :param i: the loop number for N_fill
    :param ref_file: the genome fasta file to be filled
    :param : work_dir,
    :return: the path to work dir
    """
    print(work_dir)
    tmp_dir=os.path.join(work_dir, "temp")
    tmp_ref_dir=os.path.join(tmp_dir, "ref")

    if os.path.exists(tmp_ref_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(tmp_ref_dir)

    if i==0:
        if ref_file==None:
            raise IOError("A ref file have to be provided in the first round.")
        else:
            shutil.copyfile( ref_file, (tmp_ref_dir+"/round{}.fasta".format(i)) )
    else:
        ref_last=glob((work_dir+"/*.fasta"))[-1]
        shutil.copy(ref_last, tmp_ref_dir)

    return tmp_ref_dir


@timer
def run_mapper(i,work_dir,read_list):
    """
    Considering the ref is in workdir/temp/ref
    :param work_dir:
    :param read_list:
    :param i:
    :return:
    """

    # bwa index
    ref_dir=os.path.join(work_dir, "temp/ref")
    os.chdir(ref_dir)  # set env 1st
    ref_file=myglob(ref_dir, "*.fasta")[-1]
    myexe("bwa index {0}".format(ref_file))
    print(os.listdir("."))

    # bwa mapping
    os.chdir(work_dir+"/temp")
    ## make the mapping with a large -w (to get the ungapped N) and -t (speed up)
    ## as my expreience, using a long seed and a long width will help the indel calling
    wrapper_bwamem(ref_file,read_list,prefix="round{}".format(i),core=32, w=25000, k=22)
    print(os.listdir("."))
# ----------------------------------------------------------------------------------------------------------------------
# Just test the new bam mapping data and the old

@timer
def run_nreplace(i, work_dir):
    """
    i is the loop time for the genome

    input is the ref fasta file location

    need to used 3 file:
    summary_N: N.dat
    cons: N_text.txt
    sequence_replace: replaced.fasta
    todo: should add a end point, when the N.dat do not change after a round of fix, then just stop
    """
    # dir
    ref_dir=os.path.join(work_dir, "temp/ref")

    ref=myglob(ref_dir, "*.fasta")[-1] # the latest reference should have a "larger" name
    os.chdir(work_dir)
    # summary N
    record_dict=fasta2dic(ref)
    N_roundi=summary_N(record_dict)

    # use the bam file from run_mapper
    samfile_dir=work_dir+"/temp/round{}_s.bam".format(i)

    # get N_replace.txt
    write_nreplace(record_dict= record_dict, samfile_dir=samfile_dir,
                   N_list=N_roundi, outfile=work_dir+"/N_round{}.txt".format(i),flanking=5)
    # read the file to mem
    N_replace= read_nreplace(work_dir+"/N_round{}.txt".format(i), flanking=5)

    write_nreplace_used(N_replace, outfile=work_dir+"/N_round{}_used.txt".format(i))

    if len(N_replace)<=5:
        return -1
    else:
        sequence_replace(record_dict=record_dict, N_replace=N_replace, outfile="round{}_nfill.fasta".format(i+1))
        return i+1


def run_insertion(i, workdir):
    """
    find the insertions in the long reads, fill them into the genome, several times
    :param i:
    :param workdir:
    :return:
    """
    ref_dir = os.path.join(work_dir, "temp/ref")



def clearup(work_dir):
    """
    remove the files in temp file and renew the folder
    :param work_dir:
    :return:
    """


def run_main(i=2):
    i_p=i-1


if __name__=="__main__":

    work_dir="/home/zhaolab1/myapp/LR_toolkit/test/wkdir"
    read_list=["/home/zhaolab1/myapp/LR_toolkit/test/cb12.fq"]
    #n = run_nreplace(0, work_dir=work_dir)

    for i in range(0,7):
        pre_dir_file(i,work_dir, ref_file="/home/zhaolab1/myapp/LR_toolkit/test/round0.fasta")
        run_mapper(i, work_dir,read_list )
        n=run_nreplace(i, work_dir=work_dir)
        print(n)
