#!/usr/bin/env python3

import os, gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from .BUseq import transform_seq

def msa2flat(dfasta):
    """ Transform an MSA sequence to a full protein sequence
    
    Parameter
    =========
    dfasta : dict
        a dictionary containing the MSA fasta sequences
        
    Return
    ======
    dfasta_flat : dict
        a dictionary containing the fasta sequence without MSA characters
    """
    dfasta_flat = dict()
    for prot in dfasta:
        seqobj  = dfasta[prot]
        seq = str(seqobj.seq)
        new_seq = Seq(transform_seq(seq))
        record = SeqRecord(new_seq, seqobj.id, seqobj.name, seqobj.description)
        dfasta_flat[prot] = record
    return dfasta_flat

def read_fastabio(path):
    """ read fasta using SeqIO
    """
    ext = os.path.splitext(path)[1]
    if ext in [".gzip", ".gz"]:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(path) as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return record_dict

def read_fastabioit(path):
    """ read fasta file, one sequence after an other
    """
    ext = os.path.splitext(path)[1]
    if ext in [".gzip", ".gz"]:
        with gzip.open(path, 'rt', encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record
    else:
        with open(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record

def prefix_file(path):
    return os.path.splitext(os.path.basename(path))[0]

def make_dir_ifnot(path):
    """ create a directory if not present
    """
    if not os.path.isdir(path):
        os.makedirs(path)

def read_onecol(path):
    """ read one column
    """
    data = []
    with open(path) as inf:
        for line in inf:
            data.append(line.split()[0])
    return data


def read_twocols(path):
    """ read a two column file format and store items in a
    dicionary of list
    """
    data = {}
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            data.setdefault(tmp[0], []).append(tmp[1])
    return data


def read_domains(path, **kwargs):
    """ read protein domain results
    """
    if "format" not in kwargs:
        fnct = read_pfam
    else:
        if kwargs["format"] == "pfam":
            fnct = read_pfam
        #elif kwargs["format"] == "hmmer":
            #fnct = read_hmmer
        else:
            raise ValueError("Unknown method {} passed to read_domains fucntion".format(kwargs["format"]))
    return fnct(path, kwargs)

def read_pfam(path, slct_dom=None, slct_prot=None, keep_pfamb=True, pos="ali"):
    """ read pfam domain annotation
    """
    data = {}
    with open(path) as inf:
        for line in inf:
            if line[0] != "#" and line[0] != "\n":
                tmp = line.split()
                prot = tmp[0]
                if keep_pfamb or tmp[7] != "Pfam-B":
                    if pos == "ali":
                        start, stop = int(tmp[1])-1, int(tmp[2])
                    else:
                        start, stop = int(tmp[3])-1, int(tmp[4])
                    dom = tmp[5].split(".")[0]
                    if not slct_prot or prot in slct_prot:
                        data.setdefault(prot, []).append((start, stop, dom))
    if slct_dom:
        prots = data.keys()
        for prot in prots:
            keep=False
            for start, stop, dom in data[prot]:
                if dom in slct_dom:
                    keep = True
            if not keep:
                del data[prot] 
    return data 
