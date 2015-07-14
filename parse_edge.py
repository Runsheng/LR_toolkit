#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
from utils import max_get


def edge_parse(samfile):
    """
    :param the samfile is a pysam.Alignmentfile object;
    :return: a list with (chr,start,l_or_r,seq) tuple, for the left clippling,
            the start is 0, for right clippling, the start is the length of chr
    """
    edge_list=[]
    ref_length=dict(zip(samfile.references, samfile.lengths))
    for ref in ref_length.keys():
        clipping_l=[]
        for read in samfile.fetch(ref, 0, 100):
            # cigartuples, 4 is softclipping, 5 is hardclipping,
            # TODO: the hard clipping has been ignored, can only be recovered when the fastq file directly
            if read.cigartuples[0][0]==4:
                clipping_l.append(read.seq[0:int(read.cigartuples[0][1])])
        try:
            seq_l=max_get(clipping_l) # to avoid the output of empty sequence name line
            line_l=(ref,0,"L",seq_l)
            edge_list.append(line_l)
        except Exception as e:
            print e
            print "No left clipping signal found in chr %s ." % ref

        clipping_r=[]
        for read in samfile.fetch(ref, (ref_length[ref]-100), ref_length[ref]):
            if read.cigartuples[-1][0]==4:
                clipping_r.append(read.seq[0:int(read.cigartuples[0][1])])
        try:
            seq_r=max_get(clipping_r)
            line_r=(ref, ref_length[ref]-1, "R", seq_r)
            edge_list.append(line_r)
        except Exception as e:
            print e
            print "No right clipping signal found in chr %s ." % ref

    return edge_list


if __name__=="__main__":
    samfile=pysam.AlignmentFile("./4st_split_segment/cb12nin_s.bam", "rb")
    edge_list=edge_parse(samfile)