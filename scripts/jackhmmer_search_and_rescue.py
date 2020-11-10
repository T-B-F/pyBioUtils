#!/usr/bin/env python
""" from initial sequence use jackhmmer to retrieve remote target
filter target based on the presence of the non degenerated walker A 
"""

import os, sys, argparse
import shlex, subprocess
import re
from Bio import AlignIO

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="initialseq")
    parser.add_argument("-w", action="store", dest="workdir")
    parser.add_argument("-d", action="store", dest="dbtarget")
    parser.add_argument("-n", action="store", dest="niter", type=int)
    parser.add_argument("-o", action="store", dest="output")
    params = parser.parse_args()
    return params

def run_phmmer(inpath, outpath, target):
    """ from an initial protein sequence query a database
    """
    ali_out = outpath+"_msa.fasta"
    res_out = outpath+"_res.txt"
    log_file = outpath+"_plog.txt"
    command = "phmmer -A {} -o {} {} {}".format(ali_out, res_out, inpath, target)
    #print(command)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        #try:
        subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        #except:
        #    print("Error, unable to run {}".format(command), file=sys.stderr)
        #    sys.exit(1)
    return ali_out

def run_hmmbuild(inpath, outpath):
    """ build an hmm out of a fasta file
    """
    hmmout = outpath+".hmm"
    log_file = outpath+"_hmmlog.txt"
    command = "hmmbuild {} {}".format(hmmout, inpath)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        #try:
        subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        #except:
        #    print("Error, unable to run {}".format(command), file=sys.stderr)
        #    sys.exit(1)
    return hmmout

def run_hmmsearch(inpath, outpath, target):
    """ run hmmsearch and keep aligned hits
    """
    ali_out = outpath+"_msa.fasta"
    res_out = outpath+"_res.txt"
    log_file = outpath+"_hlog.txt"
    command = "hmmsearch --noali -A {} -o {} {} {}".format(ali_out, res_out, inpath, target)
    #print(command)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        #try:
        subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        #except:
        #    print("Error, unable to run {}".format(command), file=sys.stderr)
        #    sys.exit(1)
    return ali_out

def read_results(path):
    """ read aligned sequences of hmmsearch in stockholm format.
    """
    proteins = dict()
    with open(path) as inf:
        msa = AlignIO.read(inf, format="stockholm")
    for record in msa:
        proteins[record.id] = str(record.seq)
    return proteins

def read_input(path):
    """ read initial input
    """
    proteins = set()
    with open(path) as inf:
        for line in inf:
            if line[0] == ">":
                proteins.add(line[1:].strip().split()[0])
    return proteins

def filter_results(proteins, name, workdir, mean_size):
    """ filter results based on hit coverage of initial input and on regexp presence/absence
    """
    outpath = os.path.join(workdir, name+"_filteredmsa.fasta")
    r = re.compile("G.{4}GK[TS]")#.{20,100}[YIMLFWV]{3}[YIMLFWVN]D[DE]")
    kept = set()
    with open(outpath, "w") as outf:
        for prot in proteins:
            seq_msa = proteins[prot]
            seq = seq_msa.replace("-", "").replace(".", "").upper()
            # regexp match
            match = r.search(seq)
            if not match:
                if (len(seq) / mean_size) >= 0.8:
                    outf.write(">{}\n{}\n".format(prot, seq_msa))
                    kept.add(prot.split()[0])
    return kept, outpath


def read_results_and_filter(ali_results, name, workdir, n, mean_size):
    """ apply read and filter to results
    """
    # read results
    hit_proteins = read_results(ali_results)
    # filter results based on walker A regexp
    res_proteins, filtered_ali = filter_results(hit_proteins, name+"_iter_{}".format(n), workdir, mean_size)
    return res_proteins, filtered_ali

def compute_mean_length(path):
    """ compute initial mean sequence compute_mean_length
    """
    mean_size = 0
    fasta = dict()
    with open(path) as inf:
        for line in inf:
            if line[0] == ">":
                prot = line[1:-1]
                fasta[prot] = ""
            else:
                fasta[prot] += line.strip().replace("-", "").replace(".", "")
    for prot in fasta:
        mean_size += len(fasta[prot])
    mean_size /= len(fasta)
    return mean_size


def count_differences(query_proteins, targets):
    notfound = 0
    overlapping = set()
    for prot in query_proteins:
        tmp = prot.split("/")
        name = tmp[0]
        start, stop = tmp[1].split("-")
        init_start, init_stop = int(start)-1, int(stop)
        if name in targets:
            found = False
            for new_start, new_stop in targets[name]:
                start = max(init_start, new_start)
                stop = min(init_stop, new_stop)
                diff = stop - start 
                if diff > 0:
                    c = diff / max(new_stop-new_start, init_stop-init_start)
                    if c > 0.9:
                        found = True
                        overlapping.add(prot)
            if not found:
                notfound += 1
        else:
            notfound += 1
    return notfound, overlapping

def check_set_differences(new_proteins, prev_proteins, cov=0.9):
    """ check hits, count number of new sequences and dropped sequences
    """
    init = dict()
    overlapping = set()
    for prot in prev_proteins:
        tmp = prot.split("/")
        start, stop = tmp[1].split("-")
        init.setdefault(tmp[0], list()).append((int(start)-1, int(stop)))
    
    new = dict()
    for prot in new_proteins:
        tmp = prot.split("/")
        start, stop = tmp[1].split("-")
        new.setdefault(tmp[0], list()).append((int(start)-1, int(stop)))
        
    nbnew, overlapping = count_differences(new_proteins, init)
    dropped, _= count_differences(new_proteins, init)
            
    print(nbnew, dropped, len(new_proteins), len(prev_proteins), len(overlapping))
    inboth = len(overlapping)
    if len(overlapping) == len(new_proteins):
        return True
    return False

def main():
    params = get_cmd()
   
    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)

    stop = False
    init_proteins = read_input(params.initialseq)
    init_size = len(init_proteins)
    init_length = compute_mean_length(params.initialseq)
    inpath = params.initialseq
    name = os.path.splitext(os.path.basename(inpath))[0]
    if init_size == 0:
        print("Error, no protein found in input fasta file", file=sys.stderr)
        sys.exit(0)

    elif init_size == 1:
        # single sequence
        outpath = os.path.join(params.workdir, name+"_iter_0")
        # run phmmer
        ali_results = run_phmmer(inpath, outpath, params.dbtarget)
        res_proteins, filtered_ali_results = read_results_and_filter(ali_results, name, params.workdir, 0, init_length)
        if len(res_proteins.intersection(init_proteins)) == len(res_proteins):
            # no new proteins
            stop = True
        else:
            init_proteins = res_proteins
            # output alignment of jackhmmer is on stockholm format, convert to fasta
            #with open(filtered_ali_results) as inf:
            #    msa = AlignIO.read(inf, format="stockholm")
            #with open(filtered_ali_results, "w") as ouf:
            #    AlignIO.write(msa, outf, format="fasta")
            inpath = filtered_ali_results

    if not stop:
        niter = params.niter +1 if init_size == 1 else params.niter
        for n in range(1, niter):
            outpath = os.path.join(params.workdir, name+"_iter_{}".format(n))

            # convert fasta inputfile to hmm file
            hmmpath = run_hmmbuild(inpath, outpath)

            # run hmmsearch
            ali_results = run_hmmsearch(hmmpath, outpath, params.dbtarget)
            res_proteins, filtered_ali_results = read_results_and_filter(ali_results, name, params.workdir, n, init_length)

            if check_set_differences(res_proteins, init_proteins):
                # no new proteins
                break
            if len(res_proteins) == 0:
                print("Error, no protein found", file=sys.stderr)
                break
            else:
                init_proteins = res_proteins
                # output alignment of jackhmmer is on stockholm format, convert to fasta
                #with open(filtered_ali_results) as inf:
                #    msa = AlignIO.read(inf, format="stockholm")
                #with open(filtered_ali_results, "w") as ouf:
                #    AlignIO.write(msa, outf, format="fasta")
                inpath = filtered_ali_results

    if len(init_proteins) != 0:
        if params.initialseq == inpath:
            print("No new sequences found")
        else:
            with open(params.output, "w") as outf, open(inpath) as inf:
                for line in inf:
                    outf.write(line)
    else:
        print("Exit, without proteins")
        
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()
