#!/root/anaconda3/bin/python
# -*- coding: UTF-8 -*-

## CreateTime: 2023-11-30  14:01

import os
import sys

def fiter_reads_byAS(input,as_threshold,out):
    f_in=open(input,'r')
    f_out=open(out,'w')
    lines=f_in.readlines()
    for line in lines:
        line=line.strip()
        if '@' in line:
            f_out.write(line+'\n')
        else:
            tokens=line.split('AS:i:')
            as_i=int(tokens[1].split('\t')[0])
            if as_i>as_threshold:
                f_out.write(line+'\n')
    f_in.close()
    f_out.close()
    return

def run(input,as_threshold,out):
    as_threshold=int(as_threshold)
    fiter_reads_byAS(input,as_threshold,out)