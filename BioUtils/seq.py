#!/usr/bin/env python3

msa_characters = ["-", "?", "!", "*", "."]

def transform_seq(seq):
    """ replace all characters that are not part of a protein sequence by 
    emptiness
    """
    # TODO add character checking based on ASCII code
    return "".join("" if aa in msa_characters else aa for aa in seq)
