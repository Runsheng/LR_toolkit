#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/28/16 10:22 AM
# @Author  : Runsheng
# @File    : test_mapping.py


#### define the main flow of the test run for cb4 genome
#### firstly do the mapping

import os
import shutil
import logging
import glob

from utils import myexe
from wrapper import wrapper_bwamem
from wrapper import wrapper_bam2vcf

# ----------------------------------------------------------------------------------------------------------------------
def pre_dir():
    """
    Create dirs for the mapping and working
    :return:
    """
    work_dir="./wkdir"
    tmp_dir="./wkdir/temp"
    tmp_ref_dir="./wkdir/temp/ref"

    if os.path.exists(work_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(work_dir)
        os.makedirs(tmp_dir)
        os.makedirs(tmp_ref_dir)

# ----------------------------------------------------------------------------------------------------------------------
# move the reads and the ref to wkdir
prefix=os.path.join(os.path.dirname(__file__))
index=glob.glob(prefix+"/*.fa")[0]
read_list=glob.glob(prefix+"/*.fq")

def run_mapper_caller(index,read_list):

    # set the env for the run
    prefix=os.path.join(os.path.dirname(__file__))
    work_dir=prefix+"/wkdir"
    os.chdir(work_dir)
    print(os.listdir("."))

    ## make the mapping and the vcf
    # myexe("bwa index {0}".format(index))
    # wrapper_bwamem(index,read_list)
    wrapper_bam2vcf(index,"LR_s.bam")




if __name__=="__main__":
    #pre_dir()
    run_mapper_caller(index,read_list)