#fillname: summary_N.py
#Todo, the file contians more function than expected, can be splited to utils

import cPickle as pickle
from utils import fasta2dic

def summary_N(record_dict, N_threshold=5):
    """
    Give a function to summary the "N"s in the genome, and return a bed like tuple-list to describe the "N"s,
    Output tuple format: ("chr", str,end). chr is type(str), while start and end is type(int)
    And cpickle was used to store the N_list for status count
    Note: the hard-masked genome region will be recognized as gap
    """
    N_list=[]
    for chr in record_dict.keys():
        seq=str(record_dict[chr].seq).lower()
        left=[]
        right=[]
        for i,nucl in enumerate(list(seq)):
            if nucl=="n" and seq[(i-1)]!="n":
                left.append(i)
            # the right edge can be the end of the chr
            try:
                if nucl=="n" and seq[(i+1)]!="n":
                    right.append(i)
            except Exception:
                pass
            if i==len(seq) and nucl=="n":
                right.append(i)
        # check if the left edge and right edge are paired. Normally it should be OK.
        if len(left)==len(right):
            count=0
            length_N=0
            for i in range(0,len(left)):
                N_single=(chr,left[i],right[i])
                N_list.append(N_single)
                # print some stat out
                length_N+=right[i]-left[i]+1
                if right[i]-left[i]>=N_threshold:
                    count+=1
            print "%s has %d gaps (len> %d bp), total gap length is %d (single N included)." %(chr,count,N_threshold,length_N)
        else:
            print "left", len(left), "unequal to", "right", len(right)
    # store the N_list as pickle file, leave a checkpoint:
    # with open("N_list.dat", "wb") as fp:
        #pickle.dump(N_list, fp)
    return N_list

if __name__ == '__main__':
    record_dict=fasta2dic("./4st_split_segment/cb4_insertion_filled.fasta")
    aa=summary_N(record_dict)