f=open("edge.fasta","w")

for ref in ref_length.keys():
    clipping_l=[]
    clipping_r=[]
    for read in samfile.fetch(ref, 0, 100):
        #cigartuples, 4 is softclipping, 5 is hardclipping, these sequence can extend the end of chrs
        if read.cigartuples[0][0]==4 or read.cigartuples[0][0]==5:
            clipping_l.append(read.seq[0:int(read.cigartuples[0][1])])
    for read in samfile.fetch(ref, (ref_length[ref]-100),ref_length[ref]):
        #print read.cigartuples
        if read.cigartuples[-1][0]==4 or read.cigartuples[-1][0]==5:
            clipping_r.append(read.seq[0:int(read.cigartuples[0][1])])

    try:
        max_get(clipping_l) # to avoid the output of empty sequence name line
        #print (">"+ref+"_0^0_left")
        f.write(max_get(clipping_l))
        f.write("\n")
        #print max_get(clipping_l)
        f.write(">"+ref+"_0^0_left")
        f.write("\n")
    except Exception:
        pass
        #print "No left clipping signal found."
    try:
        max_get(clipping_r)
        #print (">"+ref+"_"+str(ref_length[ref])+"^"+str(ref_length[ref])+"_right")
        #print max_get(clipping_r)
        f.write(">"+ref+"_"+str(ref_length[ref])+"^"+str(ref_length[ref])+"_right")
        f.write("\n")
        f.write(max_get(clipping_r))
        f.write("\n")
    except Exception:
        pass
        #print "No right clipping signal found."
f.close()
