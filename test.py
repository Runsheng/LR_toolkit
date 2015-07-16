#!/usr/bin/env python
# -*- coding: utf-8 -*-
from utils import myexe
import os

os.system("export LD_LIBRARY_PATH=/home/li/dna/Complete-Striped-Smith-Waterman-Library/src/")

# just test the right call
#print myexe("ls -l")

from src.ssw_wrap import Aligner
# from src.ssw_wrap import getBlastRepresentation

seq1="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGTCCCCTCTCTCT"
seq2="CGTCCCCTCTCTCT"

aligner = Aligner(seq1, report_cigar=True)
aln=aligner.align(seq2)

print aln.ref_begin,aln.ref_end,aln.query_begin,aln.query_end,aln.cigar_string

# here is a wrapper
from src.ssw_wrap import Aligner
