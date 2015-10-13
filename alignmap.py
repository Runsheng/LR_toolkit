#!/usr/bin/env python
# -*- coding: utf-8 -*-

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages

# Third part library packages

# library level import
from utils import myexe, fasta2dic

def get_edge(infile="split.fasta"):
    """
    From the stored splited segments, get the edge sequence and store them in dict
    @:param infile: a split fasta file
    @:returns
        out: dict{edge name:edge sequence}
        dict{edge_name:sequence_name}
    """
    split_dict=fasta2dic(infile)
    edge_dict={}
    name_dict={}

    for key in split_dict.keys():
        #print key
        pos_list=key.split("_")
        if "L" in key:
            len_l=int(pos_list[pos_list.index("L")+1])
            seq_l=(str(split_dict[key].seq))[:len_l]

            name_l="_".join(pos_list[0:2])+"_L_"+str(len(seq_l)) if "random" not in key else "_".join(pos_list[0:3])+"_L_"+str(len(seq_l))

            edge_dict[name_l]=seq_l
            name_dict[name_l]=key
            #print "L",name_l,len(seq_l)

        if "R" in key:
            len_r=int(pos_list[pos_list.index("R")+1])
            seq_r=str(split_dict[key].seq)[-len_r:]
            name_r="_".join(pos_list[0:2])+"_R_"+str(len(seq_r)) if "random" not in key else "_".join(pos_list[0:3])+"_R_"+str(len(seq_r))

            edge_dict[name_r]=seq_r
            name_dict[name_r]=key
            #print "R",name_r,len(seq_r)
    return edge_dict,name_dict


def run_makeblastdb():
    """
    one-time run function for the split.fasta file
    """
    out, err, proc.returncode=myexe("makeblastdb -in split.fasta -dbtype nucl")
    print out,err


def run_blast():
    """
    one-time run function tp blast the edge sequence with the split sequence
    """
    out, err, proc.returncode=myexe("blastn -query edge.fasta -db split.fasta \
                                    -num_threads 1 -max_target_seqs 1 -outfmt 6 > blastn.outfmt6")


def filter_blast(m6file="blastn.outfmt6",percent_cutoff=95,length_cutoff=100):
    """
    To filter the blast m6 table, some addtinal information need to be indentified,
    The length of the query and target sequence, the minus or plus mapping
    """
    out=[]
    with open(m6file) as f:
        lines=f.readlines()
        for line in lines:
            line_t=line.split("\t")
            name_query=line_t[0]
            name_target=line_t[1]
            percent=float(line_t[2])
            length=float(line_t[3])
            try:
                if name_query.split("_")[0]==name_target.split("_")[0] and \
                                name_query.split("_")[1]==name_target.split("_")[1]:
                    pass
                elif percent< percent_cutoff:
                    pass
                elif length < length_cutoff:
                    pass
                else:
                    out.append(line)
            except Exception as e:
                print e
                print line
    return out


def write_blast(m6file="blastn.outfmt6"):
    out=filter_blast(m6file)
    with open("blastn_filterd.m6","w") as f:
        for line in out:
            f.write(line)
    print "write blast file done!"


def parse_blast(name_dict,edge_dict,m6file="blastn_filterd.m6"):
    """
    Input: blast m6 file
    Output: the u,v digraph nodes and u->v edge
    and several status for the edge has to be generated
    for example:
        is_edge: if the mapping is in the edge of the target and/or the query
        is_unique: if the mapping only link two contig together or can link many contigs(which indicate a repeat)
                    and the two read,when to merge into consenus, has to be aligned by some other aligner so the read can overlap properly
        is_plus: check if the overlap direction is plus/plus, can repesent as (reverse(seq_name),seqname)
    """
    edge_pair=[]
    with open(m6file) as f:
        lines=f.readlines()
        for line in lines:
            line_t=line.split("\t")
            name_query=line_t[0]
            name_target=line_t[1]
            direction_query="plus" if int(line_t[7])-int(line_t[6])>0 else "minus"
            direction_target="plus" if int(line_t[9])-int(line_t[8])>0 else "minus"

            if "L" in name_query:
                edge_pair.append((name_target,name_dict[name_query],"L",direction_query,direction_target))
            elif "R" in name_query:
                edge_pair.append((name_dict[name_query],name_target,"R",direction_query,direction_target))\

    edge_f=list(set(edge_pair))
    return edge_f


def filter_edgepiar(edge_list):
    all_node=[]
    duplicate_node=[]
    # gather the duplicated nodes
    for query,target,edgetype,order_1,order_2 in edge_list:
        if query in all_node:
            duplicate_node.append(query)
        if target in all_node:
            duplicate_node.append(target)
        all_node.append(query)
        all_node.append(target)
    # gather the pair not in the duplicated nodes
    new_pair=[]
    for query,target,edgetype,order_1,order_2 in edge_list:
        if query in duplicate_node or target in duplicate_node:
            pass
        else:
            new_pair.append((query,target))
    return new_pair

def main():
    edge_dict,name_dict=get_edge(infile="./4st_split_segment/split.fasta")
    # The following code just run once
    #dic2fasta(edge_dict,out="edge.fasta")
    #print len(name_dict)
    #aa=filter_blast(m6file="blastn.outfmt6")
    #with open("blastn_filterd.m6","w") as f:
    #    for line in aa:
    #        f.write(line)
    edge_f=parse_blast(name_dict=name_dict,edge_dict=edge_dict,m6file="./4st_split_segment/blastn_filterd.m6")
    print edge_f[1]
    edge_ff=filter_edgepiar(edge_list=edge_f)
    print len(edge_ff)
    import pickle
    with open("edge_ff.dat","w") as f:
        pickle.dump(edge_ff,f)

if __name__=="__main__":
    main()