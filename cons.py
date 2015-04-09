#filename: cons.py
#Given a consesune matrix, get the cons from the sequence, used to find the replacement of the gapped "N"s

# to get per base coverage at a place
# the str of the pileup function
# Class: AlignmentFile.pileup.pileups
# Classname: AlignmentFile.pileup==pileupcolumn; AlignmentFile.pileup.pileups==pileupread
# Note: the N_list do not include the end of the gap region, so is needed to be expanded

#todo: change the marix tp pandas or nympy object
#todo: can add a function to convert the

import pysam

def con_sequence(samfile,chro,start,end):
    """
    Input: the opened pysam.Alignmentfile, and a bed interval
    Out: get the read matrix in a dict
    """
    cons_list={}
    n=0
    #take 5 mer 5' and 3' flanking sequence together with N as the input
    for pileupcolumn in samfile.pileup(chro,start-5,end+5, truncate=True):

        site_seq=""  # to avoid some UnboundLocalError in python 2.7, give an init value firstly

        for pileupread in pileupcolumn.pileups:

            if n==0: # init the first containers
                cons_list[str(pileupread.alignment.query_name)]=[]

            if pileupread.indel>0:  # for insertion
                site_seq=pileupread.alignment.query_sequence[pileupread.query_position:(pileupread.query_position+pileupread.indel+1)]
            if pileupread.is_del:  # for deletion
                site_seq="-"
            else:
                if pileupread.indel==0:  # for 1:1 matched nuclitide
                    site_seq=pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupread.indel<0:  # for the insertion nucl followed by a deletion, should be treated as the 1:1 matched
                    site_seq=pileupread.alignment.query_sequence[pileupread.query_position]

            #  need to add some new read to the dict, add more containers
            try:
                cons_list[str(pileupread.alignment.query_name)].extend(list(site_seq))
                con_length=len(cons_list[str(pileupread.alignment.query_name)])
            except KeyError:
                cons_list[str(pileupread.alignment.query_name)]=["0"]*(con_length-1)  # the current one has not been popped,so need a -1
                cons_list[str(pileupread.alignment.query_name)].extend(list(site_seq))
        n+=1  # n stand for the nucl number of the reference, not the matrix

        # Test code for the length of the cons_list
    #print chro,start,end
    return cons_list


def con_matrix(cons_list):
    """
    the samfile should be opened outside the function to avoid too much I/O
    from the same region of the bam file, get all the sequence from read for each position
    return 5 list: ATGC DEL
    """
    try:
        gaplength=len(cons_list.values()[0])

        A=[0]*gaplength
        C=[0]*gaplength
        G=[0]*gaplength
        T=[0]*gaplength
        DEL=[0]*gaplength

        for seq in cons_list.values():
            for n in range(0,len(seq)):
                if seq[n] in "Aa":
                    A[n]+=1
                if seq[n] in "Cc":
                    C[n]+=1
                if seq[n] in "Gg":
                    G[n]+=1
                if seq[n] in "Tt":
                    T[n]+=1
                if seq[n]=="-":
                    DEL[n]+=1
    # if the dict is empty, then make some empty lists directly
    except IndexError:
        gaplength=0
        A=[0]*gaplength
        C=[0]*gaplength
        G=[0]*gaplength
        T=[0]*gaplength
        DEL=[0]*gaplength

    return (DEL,A,C,G,T)

def cons(DEL,A,C,G,T):
    """
    Given: five list contains the DEL,A,C,G,T frenqency
    Return: the consensus string of the matrix
    """
    cons_string=[]
    for i in range(0,len(A)):
        col={"-":DEL[i],"A":A[i], "C":C[i], "G":G[i], "T":T[i]}
        #design model: search
        max_n=max(col.values())
        if max_n==0:  # avoid no coverage status
            pass
        else:
            for key in col.keys():
                if col[key]==max_n:
                    cons_string.append(key)
                    #just take the first as best
                    break
                    #count+=1
                #if count>1:
                    #print "Dupilicated best!"
    return "".join(cons_string)

def write_nreplace():
    # main code
    with open("N_text.txt","w") as f:
        for N_single in N_list:
            #print N_single
            chro, start,end =N_single

            matrix=con_matrix(con_sequence(samfile, chro,start,end+1))
            DEL,A,C,G,T=matrix
            sequence=cons(DEL,A,C,G,T)

            name,seq_raw=chr_select(record_dict, chro,start-5,end+1+5)

            f.write(name)
            f.write("\t")
            f.write(seq_raw)
            f.write("\t")
            f.write(sequence)
            f.write("\n")

if __name__=="__main__":
    samfile = pysam.AlignmentFile("cb12i_s.bam", "rb")
    #test code
    N_list_new=N_list[264:265]
    for N_single in N_list_new:
        print N_single
        chro, start,end =N_single
        aa=con_sequence(samfile, chro,start,end+1)
        bb=con_matrix(aa)
        DEL,A,C,G,T=bb
        sequence=cons(DEL,A,C,G,T)
        print sequence

    write_nreplace()