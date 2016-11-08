#!/usr/bin/env python
# -*- coding: utf-8 -*-

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
from collections import OrderedDict
# Third part library packages
import pysam
from utils import fasta2dic, max_get, chr_select

def parse_bps(filename):
    """
    Turn the CLC parsed breakpoint file "bp.txt" into python list
    TODO:should be modified to receive bps from different bp parser， add a para for different formats
    # need to write a parse for the CLC output of BPs
    # cal the percentage of the bps, a bp with some read covering perfectly is not a real bp.
    # the 100% single bp should be retain
    # for the pair or duplicated bps, the percentage can not be 100% (they will affect each other), however,
    # the percentage should be larger than 80% (ab cutoff)
    # generate a list containing the (start,end) tuple in order for each chr
    # this list can be directly used for the chr_select function

    Input: the output of breakpoints parsed by CLC genomic workbench 8.0 as txt file, file was sorted(as bed)
    Output: a list contains "chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage"
    Note: the bp position is 1 based
    """
    fr=open(filename,"r")
    bp_list=[]
    lines=fr.readlines()

    for n in range(1,len(lines)):  # the first line should be ignored
        linelist=lines[n].strip().split("\t")

        try:
            chro=linelist[0]
            position=linelist[1].split("^")
            start=int(position[0])
            end=int(position[1])
            l_or_r="L" if "Left" in linelist[2] else "R"

            seq=linelist[4]
            length=int(linelist[5])
            is_mapped_to_self=False if linelist[6]=="" else True
            percentage=float(linelist[-1])/(float(linelist[7])+float(linelist[8]))

        except IndexError as e:
            print e

        bp_single=chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage
        bp_list.append(bp_single)

    fr.close()
    return bp_list


def filter_bps(bp_list,percent_cutoff=0.95):
    """
    Aim: to get all the pair bps into 1 site, can treat them the same with single.
    Input: the bp_list
    Output: 3 kinds of BPs, the single, the pair and the duplication
        The "single" is the most seen one that can not be filled, usually between supercontig,
        The "duplication" can be caused by mapping, using re-alignment may change some into single,
        The "pair" is the most confidential breakpoint, usually inside the supercontig.
    """
    single=[]
    duplicate=[]
    pair=[]

    for n in range(0,len(bp_list)):
        # print n
        chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage=bp_list[n]
        try:
            chro_1,start_1,end_1,seq_1,l_or_r_1,length_1,is_mapped_to_self_1,percentage_1=bp_list[n-1] # once I made a mistake here
        except IndexError:
            print "the first bp in the full list?"
        # filter the bps: 1. some read near bps can aligned to the nearby region 2. with length less than 50
        if is_mapped_to_self or length<=50:
            pass
        else:
            if chro==chro_1 and start-start_1<=100:
                # duplication: direction the same as previous
                if l_or_r==l_or_r_1:
                    duplicate.append(bp_list[n-1])
                    duplicate.append(bp_list[n])
                else:
                    pair.append(bp_list[n-1])
                    pair.append(bp_list[n])
            # the remainings can not find a pair, and are not duplicated, marked as single,
            # some of the paired ones can also be marked, these will be retained in duplicate and pair
            else:
                if percentage>percent_cutoff:
                    single.append(bp_list[n])

    # remove the duplicated elements, the order is duplicate>pair>single
    # the intersection of these set indicate the complex duplication region
    duplicate=list(set(duplicate))
    duplicate.sort()
    pair=list(set(pair))
    pair.sort()
    single=list(set(single)-set(pair)-set(duplicate))  # almost no overlap in the 3 sets, maybe just 1 or 2, when using percentage cutoff as 0.95
    single.sort()
    return single,duplicate,pair

def filter_duplicate(duplicate,percent_cutoff=0.8):
    """
    receive a duplicate list of breakpoints, return the one in the edge,
    the format of bps is still "chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage"
    means:
    For the Left clipping (5'), chose the leftest one
    For the right clipping (3'), chose the rightest one
    todo: The pecentage should be merged
    """
    new_l=[]
    new_l.append(duplicate[0])
    for n, line in enumerate(duplicate):
        chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage=duplicate[n]
        try:
            chro_1,start_1,end_1,seq_1,l_or_r_1,length_1,is_mapped_to_self_1,percentage_1=new_l[-1]
        except IndexError:
            pass
        if chro==chro_1 and abs(start_1-start)<=100:
            if l_or_r=="L" and l_or_r_1=="L":
                percentage_n=percentage+percentage_1
                line=chro_1,start_1,end_1,seq_1,l_or_r_1,length_1,is_mapped_to_self_1,percentage_n
                new_l=new_l[:-1]
                new_l.append(line)
            elif l_or_r=="R" and l_or_r_1=="R":
                percentage_n=percentage+percentage_1
                line=chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage_n
                new_l=new_l[:-1]
                new_l.append(line)
            else:
                print "Wrong key in the l_or_r indicator"
        else:
            new_l.append(line)
    new_l_2=[]
    for line in new_l:
        if line[-1]>percent_cutoff:
            new_l_2.append(line)

    return new_l_2

def filter_pair(pair,percent_cutoff=0.8):
    """
    Receive a pair list of breakpoints, return the one in the edge,
    The pair tuple str is the same with that in single, (chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage)

    The following should be treated:
    1. some pair has dupulications and need to be merged
    2. For the 1:1 pair, (R,L) pair is OK/normal and can be treated as two single breaks, the region between the two breaks
    can be ignored as the read is rather small (< cutoff==100nt).However, remain these sequence will not harm.
    3. For the 1:1 pair, (L,R) pair will arise a new break type that can not be handled as single break, need to store
    the start and end in pair and handle them specially.
    """
    left=[]
    right=[]
    for line in pair:
        if line[4]=="L":
            left.append(line)
        else:
            right.append(line)
    left=filter_duplicate(left,percent_cutoff=0)
    right=filter_duplicate(right,percent_cutoff=0)

    print "The left break is %d and right break is %d." % (len(left),len(right))

    # the normal one is right clipping + left clipping
    pair_normal=[]
    pair_abnormal=[]
    pair_abnormal_dict={}

    if len(left)==len(right):
        print("Can be zipped together.")
        rldict=dict(zip(right,left))
        for key in sorted(rldict.keys()):
            if key[-1]+rldict[key][-1]>=percent_cutoff:
                if key[1]-rldict[key][1]<=0:
                    pair_normal.append(key)
                    pair_normal.append(rldict[key])
                else:
                    pair_abnormal.append(key)
                    pair_abnormal.append(rldict[key])
                    pair_abnormal_dict[key]=rldict[key]
    else:
        print("Can NOT be zipped together.")

    return pair_normal, pair_abnormal, pair_abnormal_dict

# import pysam
# samfile=pysam.AlignmentFile("cb12nin_s.bam", "rb")

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


def get_sequence(ref, single, duplicate_f, pair_normal, pair_abnormal, pair_abnormal_dict, edge_list):
    """
    NOTE: the bed file is 1 based! NOT 0 based.
    The str of each tuple in lists is (chro,start,end,seq,l_or_r,length,is_mapped_to_self,percentage)
    1. Select a combined (chr,start,l_or_r,seq) tuple to be used for the combination
    2. Filter all the site according to the chr, store them in a list[]
    3. Get the start position of each bp and store them in a list, then generate a tuple
        to store the starts, and get a [start,start_next] tuple for using to select chr
    """
            # the function of chr_select used the 0-based system,range as [)
            # so just use the start directly as end for previous region and
            # start for the next region from the 1 base bed file
    seq_dict=OrderedDict()
    edge_dict=OrderedDict()

    # 1. Select a combined (chr,start,l_or_r,seq) tuple to be used for the combination
    sites=[(line[0],line[1],line[4],line[3]) for line in (single+duplicate_f+pair_normal+pair_abnormal)]
    sites=sites+edge_list  # Note: do not use append, append will add one full element to form the list

    sites_pair_abnormal=[(line[0],line[1],line[4],line[3]) for line in (pair_abnormal)]
    sites_pair_abnormal_left=[(line[0],line[1],line[4],line[3]) for line in (pair_abnormal_dict.keys())]

    print "all",len(sites)
    print "carefully treated", len(sites_pair_abnormal)

    # 2. Filter all the site according to the chr, store them in a list[]
    for chro in sorted(ref.keys()):
        sites_chro=[line for line in sites if line[0]==chro]
        # hash it for easy selection
        sites_chro_dict=OrderedDict()
        for line in sites_chro:
            sites_chro_dict[line[1]]=line

        # select the start position of the breakpoint for later use
        sites_pair_abnormal_chro=[line for line in sites_pair_abnormal if line[0]==chro]
        sites_pair_abnormal_left_chro=[line for line in sites_pair_abnormal_left if line[0]==chro]

        site_start=[line[1] for line in sites_chro]
        site_start_abnormal=[line[1] for line in sites_pair_abnormal_chro]
        sites_start_abnormal_left=[line[1] for line in sites_pair_abnormal_left_chro]

        # 3. Get the start position of each bp and store them in a list, then generate a tuple as (start, next_start)
        site_start=sorted(list(set(site_start)))
        site_start_l=[0]+site_start
        site_start_r=site_start+[len(ref[chro])-1] # if use + to link two list, some order will be distruped and need to be sorted again

        site_start_l.sort()
        site_start_r.sort()
        site_start_lr=zip(site_start_l,site_start_r)
        site_start_lr.sort()

        # print site_start_lr
        for number,site_range in enumerate(site_start_lr):
            start,end=site_range

            # greedy formatting all the possibilities for the abnormal sites
            if end-start>0: # the =0 result is the normal pair and can be ignored and treated as single
                # print start,end
                if start not in site_start_abnormal and end not in site_start_abnormal:
                    name,seq=chr_select(record_dict=ref, chr=chro, start=start,end=end)
                elif end in site_start_abnormal and start not in site_start_abnormal:  # for the abnormal pair
                    end=site_start_lr[number+1][1]
                    name,seq=chr_select(record_dict=ref, chr=chro, start=start,end=end)
                elif start in site_start_abnormal and end in site_start_abnormal:
                    start=-1
                    name=""
                    seq=""
                elif start in site_start_abnormal and end not in site_start_abnormal:
                    start=site_start_lr[number-1][0]
                    name,seq=chr_select(record_dict=ref, chr=chro, start=start,end=end)

                # using the position to get the clipping information, separately
                try:
                    site_use_start = sites_chro_dict[start]
                    if site_use_start[2]=="L" and name !="":
                        name=name+"_L"+"_"+str(len(site_use_start[3]))
                        seq=seq+site_use_start[3]
                except Exception as e:
                    print e

                try:
                    site_use_end = sites_chro_dict[end]
                    if site_use_end[2]=="R" and name !="":
                        name=name+"_R"+"_"+str(len(site_use_end[3]))
                        seq=seq+site_use_end[3]
                except Exception as e:
                    print e

                print name,len(seq)
                if name !="": # to ignore the both end in the abnormal region
                    seq_dict[name]=seq
    return seq_dict


def get_edge(infile="split.fasta"):
    """
    #  Get the unalign-end from the split.fasta file
    input: a split fasta file
    out: dict{edge name:edge sequence}
         dict{edge_name:sequence_name}
    """
    split_dict=fasta2dic(infile)
    edge_dict={}
    name_dict={}

    for key in split_dict.keys():
        # print key
        pos_list=key.split("_")
        if "L" in key:
            len_l=int(pos_list[pos_list.index("L")+1])
            seq_l=(str(split_dict[key].seq))[:len_l]

            # assume that the random part of chrs is marked as "I_random"
            name_l="_".join(pos_list[0:2])+"_L_"+str(len(seq_l)) if "random" not in key else "_".join(pos_list[0:3])+"_L_"+str(len(seq_l))

            edge_dict[name_l]=seq_l
            name_dict[name_l]=key
            # print "L",name_l,len(seq_l)

        if "R" in key:
            len_r=int(pos_list[pos_list.index("R")+1])
            seq_r=str(split_dict[key].seq)[-len_r:]
            name_r="_".join(pos_list[0:2])+"_R_"+str(len(seq_r)) if "random" not in key else "_".join(pos_list[0:3])+"_R_"+str(len(seq_r))

            edge_dict[name_r]=seq_r
            name_dict[name_r]=key
            # print "R",name_r,len(seq_r)
    return edge_dict, name_dict