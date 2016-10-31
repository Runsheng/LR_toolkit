#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/28/16 10:22 AM
# @Author  : Runsheng
# @File    : test_mapping.py


#### define the main flow of the test run for cb4 genome
#### firstly do the mapping
from __future__ import print_function
import os
import shutil
import logging
import glob

# used in run_mapper
from utils import myexe
from wrapper import wrapper_bwamem
from wrapper import wrapper_bam2vcf
from utils import myglob


# only used in run_nreplace
from summary_N import summary_N
from utils import fasta2dic
from cons import write_nreplace

from sequence_replace import read_nreplace,sequence_replace


# debug functions
from utils import timer
from wrapper import show_runcmd

# ----------------------------------------------------------------------------------------------------------------------
def pre_dir():
    """
    Create dirs for the mapping and working
    :return:
    """
    root_dir=os.path.join(os.path.dirname(__file__))
    work_dir=os.path.join(root_dir, "wkdir") # remember to ignore the / in dir
    print(root_dir)
    print(work_dir)
    tmp_dir=os.path.join(root_dir, "wkdir/temp")
    tmp_ref_dir=os.path.join(root_dir, "wkdir/temp/ref")

    if os.path.exists(work_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(work_dir)
        os.makedirs(tmp_dir)
        os.makedirs(tmp_ref_dir)

    return work_dir

# ----------------------------------------------------------------------------------------------------------------------
# move the reads and the ref to wkdir
#prefix=os.path.join(os.path.dirname(__file__))
#index=glob.glob(prefix+"/*.fa")[0]
#read_list=glob.glob(prefix+"/*.fq")
#
#
#    if os.path.exists(work_dir):
#        print("Already have the dirs.")
#        pass
#    else:
#        os.makedirs(work_dir)
#        os.makedirs(tmp_dir)
#        os.makedirs(tmp_ref_dir)

# ----------------------------------------------------------------------------------------------------------------------

@show_runcmd
@timer
def run_mapper(index,read_list,i):

    # set the env for the run
    root_dir=os.path.join(os.path.dirname(__file__))
    work_dir=os.path.join(root_dir, "wkdir") # remember to ignore the / in dir
    ref_dir=os.path.join(root_dir, "wkdir/temp/ref")
    print(root_dir, work_dir, ref_dir)

    os.chdir(work_dir)
    print(os.listdir("."))

    ## make the mapping and the vcf
    # myexe("bwa index {0}".format(index))
    # myexe("mv *.fa* wkdir/temp/ref")
    # as my expreience, using a long seed and a long width will help the indel calling
    wrapper_bwamem(index,read_list,prefix="round1")

# ----------------------------------------------------------------------------------------------------------------------
# Just test the new bam mapping data and the old

@show_runcmd
@timer
def run_nreplace(i):
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
    root_dir=os.path.join(os.path.dirname(__file__))
    work_dir=os.path.join(root_dir, "wkdir") # remember to ignore the / in dir
    ref_dir=os.path.join(root_dir, "wkdir/temp/ref")
    ref=myglob(ref_dir, "*.fa")[0]

    os.chdir(work_dir)
    # summary N
    record_dict=fasta2dic(ref)
    N_round1=summary_N(record_dict)

    samfile_dir=work_dir+"/temp/round{}_s.bam".format(i)
    # get N_replace.txt
    #write_nreplace(record_dict= record_dict, samfile_dir=samfile_dir, N_list=N_round1, outfile=work_dir+"/N_round1.txt")

    N_replace= read_nreplace(work_dir+"/N_round{}.txt".format(i), flanking=10)
    if len(N_replace)<=5:
        return -1
    else:
        sequence_replace(record_dict=record_dict, N_replace=N_replace, outfile="round1_nfill.fasta")
        return i+1

if __name__=="__main__":
    #pre_dir()
    # move the reads and the ref to wkdir
    prefix = os.path.join(os.path.dirname(__file__))
    read_list = glob.glob(prefix + "/*.fq") # only used for the test
    print(read_list)

    # run_mapper_caller(index,read_list)
    run_nreplace(i=1)

