#!/usr/bin/env python
# -*- coding: utf-8 -*-
from utils import myexe
import os

#os.system("export LD_LIBRARY_PATH=/home/li/dna/Complete-Striped-Smith-Waterman-Library/ssw/")
#os.system("export LD_LIBRARY_PATH=/home/li/new/data/cb/ssw/")

# just test the right call
#print myexe("ls -l")

from ssw.ssw_wrap import Aligner
# from ssw.ssw_wrap import getBlastRepresentation

seq1="TTGAATTCCTGCCAAAAAAATTTTCTTCAAAAATTTTAATTTCCAGCCAAAATGTTTTTTTTCCGAAAATTTAAGTTTCCCGCTAAAATATTTTTCTCAAAAATTTTAATTTTCCGCCAAAAAATTTTCCAGAAAGTTTGAATTTCTTCAAAAAAT"
seq2="ATTTTCCAGAAAGTTTGAATTTCTTCAAAAAATTTGAATTCCTGCCAAAAAAATTTTCTTCAAAAATTTTAATTTCCAGCCAAAATTTTTTTTTCGAAAATTTAAATTTCCCGCTAAAATATTTTTCTCAAAAATTTTAATTTCCCGCCAAAAATATT"
aligner = Aligner(seq1, report_cigar=True)
aln=aligner.align(seq2)

print aln.ref_begin,aln.ref_end,aln.query_begin,aln.query_end,aln.cigar_string

# here is a wrapper
from ssw.ssw_wrap import Aligner
