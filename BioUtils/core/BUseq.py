#!/usr/bin/env python3

from Bio import AlignIO
import gzip, io, os

msa_characters = ["-", "?", "!", "*", "."]

def transform_seq(seq):
    """ replace all characters that are not part of a protein sequence by 
    emptiness
    """
    # TODO add character checking based on ASCII code
    return "".join("" if aa in msa_characters else aa for aa in seq)

def compute_offset_pos(seq, pos):
    """ from a sequence without gap position, computes the corresponding position in a MSA
    
    Parameters
    ==========
    seq : string
        the sequence from the MSA 
    pos : int
        the position to convert
        
    Return
    ======
    k : int
        the position in the MSA
    """
    
    nogap_seq = transform_seq(seq)
    assert(pos >= 0 and pos < len(nogap_seq))

    maps = dict()
    cnt = 0
    maxi = 0
    for i in range(len(seq)):
        if seq[i] not in msa_characters:
            maps[i-cnt] = i
            maxi = i
        else:
            cnt += 1
    return maps.get(pos, maxi)
    
    #cnt = 0
    #k = 0
    #while k<len(seq):
        #print(k, cnt, seq[k])
        #offset = 0
        #while k+offset < len(seq) and seq[k+offset] in msa_characters:
            #offset += 1
        #else:
            #cnt += 1
        #k+=offset+1
        #if cnt == pos:
            #break
    #return k
        
    #k = 0 
    #cnt = 0 if seq[k] not in msa_characters else -1
    #while cnt != pos and k < len(seq):
        #if seq[k] not in msa_characters:
            #cnt += 1
        #k += 1 
        ##print(pos, cnt, k, seq)
    #return k

def compute_revoffset_pos(seq, pos):
    """ from a MSA sequence, computes the corresponding position in the sequence without gaps
    
    Parameters
    ==========
    seq : string
        the sequence from the MSA 
    pos : int
        the position to convert
        
    Return
    ======
    k : int
        the position in the sequence without gaps
    """

    cnt = 0 
    for c in seq:
        if c in msa_characters:
            cnt += 1
    return pos - cnt

def compute_pos_seq2msa(msaseq, start, stop):
    """ from a start and stop positions of a sequence without gaps, computes the corresponding position in a MSA
    
    Parameters
    ==========
    msaseq : string
        the sequence from the MSA 
    start : int
        the position to convert
    stop : int
        the position to convert
        
    Return
    ======
    new_start : int
        the position in the MSA
    new_stop : int
        the position in the MSA
    """
    new_start = compute_offset_pos(msaseq, start)
    new_stop = compute_offset_pos(msaseq, stop-1)+1 
    return new_start, new_stop


def compute_pos_msa2seq(seq, start, stop):
    """ from a start and stop positions of MSA, computes the corresponding position in the sequence without gaps
    
    Parameters
    ==========
    seq : string
        the sequence from the MSA 
    start : int
        the position to convert
    stop : int
        the position to convert
        
    Return
    ======
    new_start : int
        the position in the sequence without gap
    new_stop : int
        the position in the sequence without gap
    """
    new_start = compute_revoffset_pos(seq, start)
    new_stop = compute_revoffset_pos(seq, stop)
    return new_start, new_stop
