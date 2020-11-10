#!/usr/bin/env python

import os, numpy as np


def profile(seq, elmt=None, window=10, stride=10):
    
    if elmt == None:
        elmt = set(seq)
    
    if isinstance(elmt, set):
        elmts = list(elmt)
    elif isinstance(elmt, list):
        elmts = elmt[:]
    elif not isinstance(elmt, list):
        elmts = [elmt]
    else:
        raise ValueError("Unknown type for parameter elmt: {}".format(elmt))
    
    # compute baseline
    mean = dict()
    size = len(seq)
    for elmt in elmts:
        cnt = 0
        for s in seq:
            if s in elmt:
                cnt += 1
        #cnt = seq.count(elmt)
        if cnt == 0:
            print("Warning, element {} not found in sequence".format(elmt))
        else:
            mean[elmt] = cnt / size
    
    # compute statistic over windows
    profiles = dict()
    for i in range(0, len(seq), stride):
        x = int(i - (window/2))
        y = int(i + (window/2))
        for elmt in mean:
            valid = 0
            cnt = 0
            for j in range(x, y):
                if j >= 0 and j < len(seq):
                    if seq[j] in elmt:
                        cnt += 1
                    valid += 1
            profiles.setdefault(elmt, list()).append((cnt/valid)-mean[elmt])
    return profiles
