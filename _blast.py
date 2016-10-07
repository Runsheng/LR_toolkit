from collections import namedtuple
blastm6_fileds = ['query', 'target', 'id', 'alignment_len', 'mistmatches', 'gap_open',
                  'q_start', 'q_end', 's_start', 's_end', "e_value", "bit_score"]
blastm6_record = namedtuple('blastm6', blastm6_fileds)


class BlastM6(object):
    """
    This BlastM6 class is used to get the many:1 hit
    from the blastp/n output with best_hit=1
    """
    def __init__(self, m6file):
        self.name=m6file

        self._blastpair_get_atr()
        self._blastpair_get()


    def _blastpair_get_atr(self):
        self.atr=[]
        with open(self.name, "r") as f:
            for line in f.readlines():
                line_t=line.strip().split("\t")
                query=line_t[0]
                blast_atr=blastm6_record(*line_t)
                self.atr.append(blast_atr)


    def _blastpair_get(self):
        """
        para:a blast m6 output file
        return: python tuples in set {(name1, name2), (name2, name3)...}
        """
        out_l=[]
        for line_t in self.atr:
            query=line_t[0]
            target=line_t[1]
            out_l.append((query, target))
        self.pairs=set(out_l)


    def filter(self, e_cut=1e-10, len_cut=50):
        """
        used to filter the blast output using atr.e_value and atr.alignment_len
        suggest: e=1e-10 and len=50 for prot, e=1e-10 and len=150 for nucl
        """
        self.atr_new=[]
        for blast_line in self.atr:
            if float(blast_line.e_value)<e_cut and float(blast_line.alignment_len)>=len_cut:
                self.atr_new.append(blast_line)

        self.atr=self.atr_new
        self._blastpair_get()