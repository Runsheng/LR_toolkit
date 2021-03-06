#!/usr/bin/env python
# -*- coding: utf-8 -*-
# filename：sequence_replace.py
# runsheng, 2015/04/09

# todo: the extend region is sovled at the split-and-link step, but if someone want to stop at this step, the extend step will be useful
# replace the sequence from the raw sequence into the new sequence
# the output is a file in lines like:
# ('I', '694217', '694228', 'GCTGGNCTATG', 'GCTGGTCTATG')==chr,start,end,seq_raw,seq_new, separated by \t

# project-level import

from utils import get_fasta_len, dic2fasta
from logger import log_detail_file
logger_f=log_detail_file("log_replace.txt")


def parse_nline(line):
    name=line.split("\t")[0]
    seq_raw=line.split("\t")[1].upper()

    # using "" as last element of empty seq_new
    try:
        seq_new=line.split("\t")[2].strip().upper()
    except IndexError:
        seq_new=""
    return name, seq_raw,seq_new


def parse_nname(name):
    chro=name.split(":")[0]
    start=name.split(":")[1].split("_")[0]
    end=name.split(":")[1].split("_")[1]
    return chro,start,end


def read_nreplace(n_text_filename, flanking=5):
    """
    read and filter the N_gap in the region, fill only the gaps jumped by the long reads
    :param n_text_filename:
    :param flanking:
    :return:
    """
    f=open(n_text_filename)

    N_replace=[]
    n_1=0
    len_1=0
    len_l_new=0

    for line in f.readlines():
        name,seq_raw,seq_new = parse_nline(line)

        # Judge the insertion credibility by 3 para: the first 5 nucls, the last 5 nucls and the length of the consensus

        # 1. "Filled", The perfect match of the sequence, the most easy way to add in
        if seq_raw[:flanking]==seq_new[:flanking] and seq_raw[-flanking:]==seq_new[-flanking:] \
                and len(seq_new)>=len(seq_raw):
            # some status
            n_1+=1
            len_1+=len(seq_raw)
            len_l_new+=len(seq_new.replace("-", ""))

            chro,start,end=parse_nname(name)
            N_single=(chro, int(start)+flanking, int(end)-flanking,
                      seq_raw[flanking:-flanking], seq_new[flanking:-flanking])
            N_replace.append(N_single)
            # some result to be known

        # 2. " Double extented or Possible double bp", The perfect match of the flanking sequence, so if there are some unalign end, the end can be used to extened
        #     the gap or used to direct the region to other place (so called bp),
        #     at this time , the seq_new is almost 10 (with some edge not well aligned, so length usually +-1)
        #  ignore at this time
        #if seq_raw[:5]==seq_new[:5] and seq_raw[-5:]==seq_new[-5:] and len(seq_new)<len(seq_raw):
        #    if len(seq_new) <14:  #all length is less than 14
                #print name,seq_raw,seq_new,len(seq_raw),len(seq_new)

        # 3. "Single extended" or "single breakpoint"
            # ignored at this time

    print("In total, %d gaps with %d bps were recorded to be filled." % (n_1, len_1))
    print("The sequences after fill will be %d long" % len_l_new)
    f.close()
    return N_replace


def write_nreplace_used(N_replace, outfile):
    """
    write a small list for used N_replace file in this round's fix
    :param N_replace:
    :return: generate a file
    """
    fw=open(outfile, "w")
    for line in N_replace:
        for i,item in enumerate(line):
            if i==len(line)-1:
                fw.write(str(item))
                fw.write("\n")
            else:
                fw.write(str(item))
                fw.write("\t")

    fw.close()


def sequence_replace(record_dict, N_replace, outfile="replaced.fasta"):
    """
    @para refdic: the reference name and sequence in dict format
    @para N_test_filename, the "N_text.txt" file, generated by the previous function
    @para outfile, the filename to be written with the replaced fasta file
    """
    new_dict={}

    for name in record_dict.keys():
        # the orign sequence, all in upper
        seq_chro=str(record_dict[name].seq).upper()

        # get fill inormation for one chr
        subreplace={}
        for N_single in N_replace:
            chro, start, end, seq_raw, seq_new=N_single
            if name==chro:
                subreplace[(int(start),int(end))]=N_single
        # get the start and end, change them to int and sort
        cutsite=[]
        for key in subreplace.keys():
            chro, start, end, seq_raw, seq_new=subreplace[key]
            cutsite.append((int(start),int(end)))
        cutsite.sort()
        # todo: change to logger
        # print("Total %d gaps filled for this chro" % len(cutsite))

        seq_chro_new=[]

        # for a chr without N sites
        if len(cutsite)==0:
            seq_chro_new_str=seq_chro
        else:
            for i in range(0,len(cutsite)):
                # start
                if i==0:
                    i_start,i_end=cutsite[i]
                    seq_1=seq_chro[:i_start]

                    seq_2_chro=seq_chro[i_start:i_end]
                    seq_2_raw=subreplace[(i_start,i_end)][3]

                    seq_2_new=subreplace[i_start,i_end][4]

                    if seq_2_chro==seq_2_raw:
                        seq_chro_new.append(seq_1)
                        seq_chro_new.append(seq_2_new)
                    else:
                        print("Unequal length of stored gap and actual gap position, check the reference sequence!")

                # common
                if 0<i<(len(cutsite)-1):
                    i_start,i_end=cutsite[i-1]
                    i2_start,i2_end=cutsite[i]

                    seq_1=seq_chro[i_end:i2_start]

                    seq_2_chro=seq_chro[i2_start:i2_end]
                    seq_2_raw=subreplace[(i2_start,i2_end)][3]

                    seq_2_new=subreplace[i2_start,i2_end][4]

                    if seq_2_chro==seq_2_raw:
                        seq_chro_new.append(seq_1)
                        seq_chro_new.append(seq_2_new)
                    else:
                        print("Unequal length of stored gap and actual gap position, check the reference sequence!")
                # end
                if i==len(cutsite)-1 and i!=0:
                    i_start,i_end=cutsite[i-1]
                    i2_start,i2_end=cutsite[i]

                    seq_1=seq_chro[i_end:i2_start]

                    seq_2_chro=seq_chro[i2_start:i2_end]
                    seq_2_raw=subreplace[(i2_start,i2_end)][3]
                    seq_2_new = subreplace[i2_start, i2_end][4]
                    seq_3=seq_chro[i2_end:]


                    if seq_2_chro == seq_2_raw:
                        seq_chro_new.append(seq_1)
                        seq_chro_new.append(seq_2_new)
                        seq_chro_new.append(seq_3)
                    else:
                        print("Unequal length of stored gap and actual gap position, check the reference sequence!")

                if i==len(cutsite)-1 and i==0:
                    i2_start,i2_end=cutsite[i]

                    seq_2_chro=seq_chro[i2_start:i2_end]
                    seq_2_raw=subreplace[(i2_start,i2_end)][3]
                    seq_2_new = subreplace[i2_start, i2_end][4]
                    seq_3=seq_chro[i2_end:]


                    if seq_2_chro == seq_2_raw:
                        seq_chro_new.append(seq_3)
                    else:
                        print("Unequal length of stored gap and actual gap position, check the reference sequence!")
            # the deletions is still in "-", remove them
            seq_chro_new_str="".join(seq_chro_new).replace("-","")

        # todo: change to logger
        logger_f.debug("Length after fill: %d; length before fill: %d." % (len(seq_chro_new_str), len(seq_chro)))

        new_dict[name]=seq_chro_new_str
        # write the sequence

    dic2fasta(new_dict, outfile)



def check_post_replace(fasta_before, fasta_after, delta):

    lb=get_fasta_len(fasta_before)
    la=get_fasta_len(fasta_after)

    if la==lb+delta:
        return True
    return False



if __name__=="__main__":
    # test code
    N_replace=read_nreplace(n_text_filename="/home/zhaolab1/myapp/LR_toolkit/test/wkdir/N_round5.txt")
    #print N_replace[2]

    # main code
    #sequence_replace()
