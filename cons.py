# to get per base coverage at a place
# the str of the pileup function
# Class: AlignmentFile.pileup.pileups
# Classname: AlignmentFile.pileup==pileupcolumn; AlignmentFile.pileup.pileups==pileupread
# Note: the N_list do not include the end of the gap region, so is needed to be expanded

#filename: cons.py
#Given a consesune matrix, get the cons from the sequence, used to find the replacement of the gapped "N"s

#todo: change the marix tp pandas or nympy object
#todo: add the detection of insertion

import pysam

def con_sequence(samfile,chro,start,end):
    """
    from the pysame.Alignmentfile, get the reads in the N region, include the insertion and deletion
    """
    cons_list={}
    n=0
    #take the 5 mer 5' and 3' flanking sequence as the input
    for pileupcolumn in samfile.pileup(chro,start-5,end+5, truncate=True):
        for pileupread in pileupcolumn.pileups:
            if n==0: #init all the containers
                cons_list[str(pileupread.alignment.query_name)]=[]

            if pileupread.indel>0:  # for insertion
                site_seq=pileupread.alignment.query_sequence[pileupread.query_position:(pileupread.query_position+pileupread.indel+1)]
            if pileupread.is_del:  #for deletion
                site_seq="-"
            else:
                if pileupread.indel==0:  # for 1:1 match N:N
                    site_seq=pileupread.alignment.query_sequence[pileupread.query_position]

            cons_list[str(pileupread.alignment.query_name)].extend(list(site_seq))
        n+=1
    #print n
    return cons_list


def con_matrix(samfile,chro,start,end):
    """
    the samfile should be opened outside the function to avoid too much I/O
    from the same region of the bam file, get all the sequence from read for each position
    return 5 list: ATGC DEL
    """
    #initiate the matrix，need to be reasigned
    gaplength=end-start # note the end need to be +1 in the input, from the N_list
    A=[0]*gaplength
    C=[0]*gaplength
    G=[0]*gaplength
    T=[0]*gaplength
    DEL=[0]*gaplength

    #write values to the matrix
    #直接选头尾所在的序列位置？
    n=0
    for pileupcolumn in samfile.pileup(chro,start,end, truncate=True):
        for pileupread in pileupcolumn.pileups:
            site_seq=pileupread.alignment.query_sequence[pileupread.query_position]
            #print pileupread.indel
            #print pileupread.query_position
            #test if it is del first
            if pileupread.is_del:
                DEL[n]+=1
            elif site_seq in "Aa":
                A[n]+=1
            elif site_seq in "Cc":
                C[n]+=1
            elif site_seq in "Gg":
                G[n]+=1
            elif site_seq in "Tt":
                T[n]+=1

        coverage=pileupcolumn.n  # for each site, get the coverage, the coverage should be expained fullly by the deletion and the nuclitides
        #add a test of coverage
        if coverage==(DEL[n]+A[n]+C[n]+G[n]+T[n]):
            pass
        else:
            print "Coverage count error at %s:%s" % (chro,pileupcolumn.pos)
        n+=1

    return (DEL,A,C,G,T)

def cons(DEL,A,C,G,T):
    """
    Given: four list contains the DEL,A,C,G,T frenqency
    Return: the consensus string of the matrix
    """
    cons_string=[]
    for i in range(0,len(A)):
        col={"-":DEL[i],"A":A[i], "C":C[i], "G":G[i], "T":T[i]}
        #design model: search
        max_n=max(col.values())
        if max_n==0:  #avoid no coverage status
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

if __name__=="__main__":
    samfile = pysam.AlignmentFile("cb12i_s.bam", "rb")
    N_list_new=N_list[136:137]
    for N_single in N_list_new:
        #print N_single
        chro, start,end =N_single
        aa=con_sequence(samfile, chro,start,end+1)
        DEL,A,C,G,T=con_matrix(samfile, chro,start,end+1) #note the +1 is still needed
        sequence=cons(DEL,A,C,G,T)
        print sequence