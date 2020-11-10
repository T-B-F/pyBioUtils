#!/usr/bin/env python
""" cut a large fasta file into multiple smaller files
"""

import os, sys, argparse
from myPython.my_utils import prefix_file, read_fastabioit

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="pathin")
    parser.add_argument("-o", action="store", dest="pathout")
    parser.add_argument("-s", action="store", dest="size", type=int)
    parser.add_argument("-n", action="store", dest="nbfile", type=int)
    params = parser.parse_args()
    return params

def split_filelist(pathin, pathout, size):
    name = prefix_file(pathin)
    cnt_file = 0
    cnt = 0
    with open(pathin) as inf:
        outfile = open(os.path.join(pathout, name+".{}.dat".format(cnt_file)), "w")
        for line in inf:
            outfile.write(line)
            cnt += 1 
            if cnt == size:
                outfile.close()
                cnt_file += 1
                cnt = 0
                outfile = open(os.path.join(pathout, name+".{}.dat".format(cnt_file)), "w")
        outfile.close()

def main():
    
    params = get_cmd()
    if not os.path.isdir(params.pathout):
        os.makedirs(params.pathout)

    if params.size:
        split_filelist(params.pathin, params.pathout, params.size)
    elif params.nbfile:
        with open(params.pathin) as inf:
            nb_line = sum(1 for line in inf)
        size = nb_line // params.nbfile
        if nb_line % params.nbfile != 0:
            size += 1
        split_filelist(params.pathin, params.pathout, size)

    sys.exit(0)
        

if __name__ == "__main__":
    main()
