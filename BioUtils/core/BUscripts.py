#!usr/bin/env python

import BioUtils
from Bio import SeqIO
import numpy as np
import os, sys, argparse

def get_cmd_split():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input fasta file")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output fasta file")
    parser.add_argument("-p", "--percentage", action="store", dest="percentage", help="percentage of proteins", type=float)
    parser.add_argument("-n", "--number", action="store", dest="number", help="first n proteins", type=int)
    params = parser.parse_args()

    if params.percentage is not None and params.number is not None:
        raise ValueError("Either percentage or number argument must be set, not both")
    if params.percentage is None and params.number is None:
        raise ValueError("Either percentage or number argument must be set, not both")

    return params

def fasta_split():
    params = get_cmd_split()

    if params.percentage:
        dfasta = BioUtils.BUio.read_fastabio(params.input)
        proteins = list(dfasta.keys())
        length = len(proteins)
        nb = int((length * params.percentage) / 100)
    else:
        nb = params.number

    with open(params.output, "w") as outf:
        cnt = 0
        for rec in BioUtils.BUio.read_fastabioit(params.input):
            if cnt < nb:
                SeqIO.write(rec, outf, "fasta")
                cnt += 1
    
    sys.exit(0)

###### shuffle

def get_cmd_shuffle():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input fasta file")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output fasta file")
    params = parser.parse_args()
    return params

def fasta_shuffle():
    params = get_cmd_shuffle()
    
    dfasta = BioUtils.BUio.read_fastabio(params.input)
    proteins = list(dfasta.keys())
    np.random.shuffle(proteins)
    with open(params.output, "w") as outf:
        for prot in proteins:
            SeqIO.write(dfasta[prot], outf, "fasta")
    
    sys.exit(0)

###### convert

def get_cmd_convert():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input sequence file")
    parser.add_argument("-o", "--outdir", action="store", dest="output", help="output sequence file")
    parser.add_argument("--input_format", action="store", dest="input_format", help="input format")
    parser.add_argument("--output_format", action="store", dest="output_format", help="output format ")
    params = parser.parse_args()

    return params

def fasta_convert():
    params = get_cmd_convert()
    
    with open(params.input, "r") as inf, open(params.output, "w") as outf:
        for rec in SeqIO.parse(inf, params.input_format):
            SeqIO.write(rec, outf, params.output_format)
    
    sys.exit(0)

###### batch

def get_cmd_batch():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input fasta file")
    parser.add_argument("-o", "--outdir", action="store", dest="output", help="output directory")
    parser.add_argument("-p", "--output", action="store", dest="prefix", help="prefix of output file")
    parser.add_argument("-b", "--batch_size", action="store", dest="batch_size", help="size of batch", type=int)
    parser.add_argument("-n", "--nb_batch", action="store", dest="nb_batch", help="number of batch", type=int)
    params = parser.parse_args()

    if params.batch_size is not None and params.nb_batch is not None:
        raise ValueError("Either batch_size or nb_batch argument must be set, not both")
    if params.batch_size is None and params.nb_batch is None:
        raise ValueError("Either batch_size or nb_batch argument must be set, not both")

    return params

def fasta_batch():
    params = get_cmd_batch()
    
    if not os.path.isdir(params.output):
        os.makedirs(params.output)

    if params.percentage:
        dfasta = BioUtils.BUio.read_fastabio(params.input)
        proteins = list(dfasta.keys())
        length = len(proteins)
        nb = int((length * params.nb_batch) / 100)
    else:
        nb = params.number

    idx = 1
    outpath = os.path.join(params.output, params.prefix+"_{}.fasta".format(idx))
    outf = open(outpath, "w")
    cnt = 0
    for rec in BioUtils.BUio.read_fastabioit(params.input):
        if cnt >= nb:
            outf.close()
            idx += 1
            outpath = os.path.join(params.output, params.prefix+"_{}.fasta".format(idx))
            outf = open(outpath, "w")
            cnt = 0
        SeqIO.write(rec, outf, "fasta")
        cnt += 1
    
    sys.exit(0)