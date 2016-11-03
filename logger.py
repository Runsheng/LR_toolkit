#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/6/16 12:47 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : logger.py


import logging


def log_summary():
    """
    # define only the summary logger first
    """
    logger = logging.getLogger('summary')
    logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def log_detail_file(filename):
    """
    # define only the summary logger first
    """
    logger = logging.getLogger('details')
    logger.setLevel(logging.DEBUG)

    ch = logging.FileHandler(filename=filename)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def test_summary():
    logger=log_summary()
    logger.info("This is just a test, print number 100 as float"+","+str(float(1)*100))
    logger.info("In total %d line, %s" % (2, "just a test number"))


def test_detail_file():
    logger=log_detail_file("test.txt")
    logger.debug("This is just a test, print number 100 as float"+","+str(float(1)*100))
    logger.debug("In total %d line, %s" % (2, "just a test number"))


if __name__=="__main__":
    import os
    os.chdir("./test")
    test_summary()
    test_detail_file()
