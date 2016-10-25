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

from utils import myexe
from wrapper import wrapper_bwamem
from wrapper import wrapper_bam2vcf
from utils import myglob

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


def run_mapper_caller(index,read_list):

    # set the env for the run
    root_dir=os.path.join(os.path.dirname(__file__))
    work_dir=os.path.join(root_dir, "wkdir") # remember to ignore the / in dir
    print(root_dir)
    print(work_dir)

    os.chdir(work_dir)
    print(os.listdir("."))

    ## make the mapping and the vcf
    # myexe("bwa index {0}".format(index))
    wrapper_bwamem(index,read_list,prefix="round1_") # as my expreience, using a long seed and a long width will help the indel calling
    #wrapper_bam2vcf(index,"LR_s.bam")

# ----------------------------------------------------------------------------------------------------------------------
# Just test the new bam mapping data and the old

def run_breplace():
    """
    need to used 3 file:
    summary_N: N.dat
    cons: N_text.txt
    sequence_replace: replaced.fasta
    todo: should add a end point, when the N.dat do not change after a round of fix, then just stop
    """

    record_dict=fasta2dic("./test/")
    aa=summary_N(record_dict)
    myglob()



if __name__=="__main__":
    #pre_dir()
    # move the reads and the ref to wkdir
    prefix = os.path.join(os.path.dirname(__file__))
    index = glob.glob(prefix + "/*.fa")[0]
    read_list = glob.glob(prefix + "/*.fq")
    print(read_list, index)
    run_mapper_caller(index,read_list)

