#!usr/bin/env python

import BioUtils
from Bio import SeqIO, AlignIO, Alphabet
import numpy as np
import os, sys, argparse, gzip, tempfile

###### oneline

def get_cmd_oneliner():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input fasta file")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output fasta file")
    params = parser.parse_args()

    return params

def fasta_oneline():
    params = get_cmd_oneliner()

    if params.input != params.output:
        with open(params.output, "w") as outf:
            for rec in BioUtils.BUio.read_fastabioit(params.input):
                outf.write(">{} {}\n{}\n".format(rec.id, rec.description, str(rec.seq)))
    else:
        tmp, tmppath = tempfile.mkstemp()
        with os.fdopen(tmp, 'w') as outf:
            for rec in BioUtils.BUio.read_fastabioit(params.input):
                outf.write(">{} {}\n{}\n".format(rec.id, rec.description, str(rec.seq)))
                
        with open(tmppath, "r") as inf, open(params.output, "w") as outf:
            for line in inf:
                outf.write(line)
    
    sys.exit(0)



###### split

def get_cmd_split():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input fasta file")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output fasta file")
    parser.add_argument("-p", "--percentage", action="store", dest="percentage", help="percentage of proteins", type=float)
    parser.add_argument("-n", "--number", action="store", dest="number", help="first n proteins", type=int)
    parser.add_argument("-s", "--shuffle", action="store_true", dest="shuffling", help="shuffle sequences before split")
    params = parser.parse_args()

    if params.percentage is not None and params.number is not None:
        raise ValueError("Either percentage or number argument must be set, not both")
    if params.percentage is None and params.number is None:
        raise ValueError("Either percentage or number argument must be set, not both")

    return params

def fasta_split():
    params = get_cmd_split()

    if params.shuffling:
        proteins = list()
        for rec in BioUtils.BUio.read_fastabioit(params.input):
            proteins.append(rec.id)
        np.random.shuffle(proteins)
        if params.percentage:
            length = len(proteins)
            nb = int((length * params.percentage) / 100)
        else:
            nb = params.number
            
        proteins = set(proteins[:nb])
        with open(params.output, "w") as outf:
            for rec in BioUtils.BUio.read_fastabioit(params.input):
                if rec.id in proteins:
                    SeqIO.write(rec, outf, "fasta")
                
    else:
        if params.percentage:
            length = 0
            ext = os.path.splitext(params.input)
            if ext in [".gzip", ".gz"]:
                with gzip.open(params.input, "rt") as inf:
                    for line in inf:
                        if line[0] == ">":
                            length += 1
            else:
                with open(params.input, "r") as inf:
                    for line in inf:
                        if line[0] == ">":
                            length += 1
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
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", help="prefix of output file")
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
    
    
###### msa convert

def get_cmd_msaconvert():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input sequence file")
    parser.add_argument("-o", "--outdir", action="store", dest="output", help="output sequence file")
    parser.add_argument("-a", "--alphabet", action="store", dest="alphabet", help="sequences are either dna, rna or protein sequences")
    parser.add_argument("--input_format", action="store", dest="input_format", help="input format")
    parser.add_argument("--output_format", action="store", dest="output_format", help="output format ")
    params = parser.parse_args()

    return params

def msa_convert():
    params = get_cmd_msaconvert()
    
    alphabet = {"dna": Alphabet.generic_dna,
                "rna": Alphabet.generic_rna,
                "protein": Alphabet.generic_protein}
            
    with open(params.input, "r") as inf, open(params.output, "w") as outf:
        if params.alphabet is not None:
            msa = AlignIO.read(inf, params.input_format, alphabet = alphabet[params.alphabet])
        else:
            msa = AlignIO.read(inf, params.input_format)
        AlignIO.write(msa, outf, params.output_format)
    
    sys.exit(0)

###### clean msa col

def get_cmd_cleancol():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", dest="input", help="input MSA file")
    parser.add_argument("--input_format", action="store", dest="input_format", help="MSA format", default="fasta")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output file")
    params = parser.parse_args()

    return params

def msa_cleancol():
    params = get_cmd_cleancol()
    
    with open(params.input, "r") as inf:
        msa = AlignIO.read(inf, params.input_format)
    
    nb_char = len(msa[0])
    nb_seq = len(msa)
    
    to_remove = list()
    for i in range(nb_char):
        col = msa[:, i]
        nb_gap = col.count("-") + col.count(".")
        if nb_gap == nb_seq:
            to_remove.append(i)
    
    for offset, i in enumerate(to_remove):
        msa = msa[:, :i-offset] + msa[:, i+1-offset:]
          
    with open(params.output, "w") as outf:
        AlignIO.write(msa, outf, params.input_format)
    
    sys.exit(0)

