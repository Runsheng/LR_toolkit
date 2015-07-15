#!/usr/bin/env python
# -*- coding: utf-8 -*-
from utils import myexe
import os

os.system("export LD_LIBRARY_PATH=/home/li/dna/Complete-Striped-Smith-Waterman-Library/src/")

# just test the right call
#print myexe("ls -l")

from src.ssw_wrap import Aligner
from src.ssw_wrap import getBlastRepresentation

seq1="AATCGTCCCCTTTTT"
seq2="AGTCCCCTTTTT"

aligner = Aligner(seq1, report_cigar=True, report_secondary=True, match=4,mismatch=-1,gap_open=-1)
aa=aligner.align(seq2)



print getBlastRepresentation(seq1,seq2, aa.cigar_string)

