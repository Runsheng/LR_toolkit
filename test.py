#! /usr/bin/env python
from utils import myexe

print myexe("ls -l")

from src.ssw_wrap import Aligner
from src.ssw_wrap import getBlastRepresentation
aligner = Aligner("AGTCGT", report_cigar=True, report_secondary=True)
aa=aligner.align("AGTC")
print aa.ref_begin,aa.ref_end,aa.query_begin,aa.query_end,aa.cigar_string

print getBlastRepresentation(aa)