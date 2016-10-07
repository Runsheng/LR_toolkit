#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/6/16 11:37 AM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : gff.py

# sys level import
from collections import namedtuple, OrderedDict
import os

# self import
from logger import log_summary

# set logger
import logging
logger=log_summary()


# set nametuple for gff class
GFF_FIELD = ['seqid', 'source', 'type', 'start', 'end',
                 'score', 'strand', 'frame', 'attributes']

GFF_record = namedtuple("gff", GFF_FIELD)


class GFF(object):
    """
    A general class to parse the gff file
    For speed issue, just read all lines into mem
    Todo: re-write as generator, reduce the flexibility but improve the mem use
    """

    def __init__(self, gfffile):
        self.filename=gfffile

        self._gff_string = []
        self.gff_to_list() # here use list for the flexibility instead of generator

        self.gff_records=[] # hold multiple GFF_record in a list
        self.parse_gff_line()

    def gff_to_list(self):
        """
        :return: gfflines list
        """
        with open(self.filename, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != "#" and len(line) > 18: # ignore the # and small lines
                    self._gff_string.append(line)

    def parse_gff_line(self):
        """
        :return: gff_records list, contains 9 col nametuple
        """

        for line in self._gff_string:
            line_l=line.strip().split("\t")

        # test if the line is 9 cols
            if len(line_l) != len(GFF_FIELD):
                logging.error(line)
                logging.error("{0} != {1}".format(len(line_l), len(GFF_FIELD)))
                pass

        # in case some gff have quotes in the cols, replace them
            else:
                data = {
                    'seqid': None if line_l[0] == '.' else line_l[0].replace('"', ''),
                    'source': None if line_l[1] == '.' else line_l[1].replace('"', ''),
                    'type': None if line_l[2] == '.' else line_l[2].replace('"', ''),
                    'start': None if line_l[3] == '.' else int(line_l[3]),
                    'end': None if line_l[4] == '.' else int(line_l[4]),
                    'score': None if line_l[5] == '.' else float(line_l[5]),
                    'strand': None if line_l[6] == '.' else line_l[6].replace('"', ''),
                    'frame': None if line_l[7] == '.' else line_l[7].replace('"', ''),
                    'attributes':self._parse_attributes(line_l[8])
                }

                self.gff_records.append(GFF_record(**data))

    @staticmethod
    def _parse_attributes(attr_string):
        """
        used to parse the GFF attributes field (col 9) into a dict
        example: ID=ctg7180000006767:hsp:104246:3.10.0.0;Parent=ctg7180000006767:hit:31434:3.10.0.0; \
        Target=CBG07519 47 108;Length=292;Gap=F2 M62 R1
        return: [('ID', 'ctg7180000006767:hsp:104246:3.10.0.0'), ('Parent', 'ctg7180000006767:hit:31434:3.10.0.0'), \
         ('Target', 'CBG07519 47 108'), ('Length', '292'), ('Gap', 'F2 M62 R1')]))
        """
        attr_dic=OrderedDict()

        if attr_string == ".":
            return attr_dic
        else:
            attr_dic = OrderedDict()
            for attribute in attr_string.strip().split(";"):
                if len(attribute)>0:
                    elems = attribute.strip().split('=')
                    key = elems[0]
                    value = ' '.join(elems[1:])
                    if value[0] == '"':
                        value = value[1:]
                    if value[-1] == '"':
                        value = value[0:-1]
                    attr_dic[key] = value
        return attr_dic

    @staticmethod
    def split(self, chro, splitpoint):
        # todo: leave the api to join
        """
        :param: chro, define the chro to be split in the
        :param: splitpoint, a list contains the 1-based points to break, for example, [600, 1200, 1900]
        :return: A dict, {name1:GFF, name2:GFF}, name use the raw seqID and the
        """
        for record_one in self.gff_records:
            if

    def join(self, nextgff):
        # todo: make this used as a suffix
        """
        :param nextgff: another GFF instance
        :return:
        """
        pass



class GFFSingle(object):
    pass


def test():
    g1=GFF("./test/example/ctg7180000006767.gff")
    print g1
    for i in g1.gff_records:
        print i

if __name__=="__main__":
    test()
    print "done"
