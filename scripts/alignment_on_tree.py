#!/usr/bin/env python
""" automatically build a taxonomically tree on protein sequence
"""

import os, sys, argparse
import subprocess, shlex
import ete3
import time
from ete3 import NCBITaxa
from ete3 import SeqMotifFace, TreeStyle, NodeStyle, add_face_to_node

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfile")
    parser.add_argument("--prot2taxid", action="store", dest="knownprot2taxid")
    parser.add_argument("--oali", action="store", dest="alitree")
    parser.add_argument("--otax", action="store", dest="sptree")
    parser.add_argument("-t", action="store", dest="outtaxid")
    params = parser.parse_args()
    return params


def read_inputfile(path):
    fasta = ""
    proteins = dict()
    with open(path) as inf:
        for line in inf:
            fasta += line.upper()
            if line[0] == ">":
                prot = line[1:].strip().split("/")[0]
                proteins.setdefault(prot, list()).append(line[1:].strip())
    return proteins, fasta

def parse_results(res):
    lines = res.split("\n")
    taxid = None
    spname = None
    for line in lines:
        line = line.strip()
        if line.startswith("<TSeq_taxid>"):
            taxid = line.lstrip("<TSeq_taxid>").rstrip("</TSeq_taxid>")
        if line.startswith("<TSeq_orgname>"):
            spname = line.lstrip("<TSeq_orgname>").rstrip("</TSeq_orgname>")
    if taxid is None or spname is None:
        raise ValueError("Unable to find taxonomic or sp name id in:\n\n{}".format(res))
    return taxid, spname
            
    
def retrieve_taxid(proteins, taxid2prot, taxid2sp, prot2taxid):
    """
    (pydev3) tristanbf@lcqb0011:~/Work/postdoc_impmc/SAMD9L$ curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=XP_017458242&amp;rettype=fasta&amp;retmode=xml"
    <?xml version="1.0" ?>
    <!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">
    <TSeqSet>
    <TSeq>
    <TSeq_seqtype value="protein"/>
    <TSeq_gi>1046831124</TSeq_gi>
    <TSeq_accver>XP_017458242.1</TSeq_accver>
    <TSeq_taxid>10116</TSeq_taxid>
    <TSeq_orgname>Rattus norvegicus</TSeq_orgname>
    <TSeq_defline>PREDICTED: LOW QUALITY PROTEIN: sterile alpha motif domain-containing protein 9-like [Rattus norvegicus]</TSeq_defline>
    <TSeq_length>1562</TSeq_length>
    <TSeq_sequence>MDRHVTQPKLIKDWTKEHVRKWITEDLKIDEKYAQIVFDEEVTGMVLQELTEKDLREMGLPRGPALLIKHMYNKLISSPEGHNQDSRQLNNKTLSINEQPKKSNSEEENSISSNSDLGLRETGQNEERGTSLLKENTLGDVETKDMKGNKPKSEQMSCMPHPFNFAHDANRYIEHSILRVAETGPLNLIDPVHEFKAFTNTETAQEQMKMKFSNETFRFAAACMNTRTNGTIHFGVKDKPHGEIVGVQVTSKDIFINHFNVMITKYFEDSDIHEARACIREPRFVEVLLLNNTQSNRFVIEVDVIPKHSICQEKYFYVMLQTCTGTTWKQSKDTSLFVREGASSRDILGNPKQRDREFKKFLENLKMSIASRKAAEECMVVSKKDSEGLKLSKLLTRHQGSLDKSYYDWYILVTNSCVPTQMEHLDFIKEIKLFAVLDFDPYSHIQGVFKAYRESRVANLHLPTHYEETTTIEEKISTLKLYEQPSWIFCNGRVDLSCQPLEPHLWQRDRASGVRKLISFLTGGNIIERGKVLVVFLLLSPVENQKDPLLETFCAFYQVFNGMDNMLCICVNPSIYQQWSDLLQVRLEMKDDLAKHSISTLNIELVNSTILKLKSVIQSSRRFLPSCGSSSVILEKMDEDIMSALEILCENECRDTDIEKDESQFLEFKKLREEHFYRGGRVSWWNFYFSSENYSSAFVKRDNFEELTTLIQQCADSPNPVFAKFINMYHHPGCGGTTLAMHVLWDLKQKFRCAVLKNKATDFIEIGEQVSKLISYKASSHQDYIPVLLLVDDFEEQEDTYILQSAINSFIAHKGVRYEKTLVIILNCMRSQNPDESAKLADSISLKYQLSPKEKRAFDAKLQEIEKEYKDCENFYSFMILKDNFDTTYIKNVVKNTLKDLDAHSRKAQLISYLALLNTFVIDSTISVSQCEIFLGITYRNKRGKLETVEDNMGTYSTLLIRTEVSDYGKYAGIRITHPLIAIHCLKELEMKYGMDRCHIALNMLEENVFYNSGIGRDKFKHDVQTLLLTRQCKEHGGEIGTLFSPFMEELQDEETEKVLTAGSNRFPQNTFIRRALARHFYLKEKNFSTALVWANQAKRKAPMNLYISDKLGRVYKRDIRCWLGKNNTYRSLSVDDLTRFLDVGEKASKAFKESQEQTDSKDYQTEFWSPQMFRRRTLNTAGFFGEIEVGLDTIQLRQLTPPFHKENEMSKESMAQFSSGKGTIPPEPKVLSNFITFLQSLKRHFDFFLDYMGLLKPRNTPKELTELSLSKKVSRCFKKYAERFCQELQGKEDLLLQEENCRKRIKAWRADTFSGLLEYLNPNHKEANNMENIVEHYTFLLQRTLNKQLSKALIKDTQNFILANIIISCLKPSSKHILPFSTLKKKLREVLQIVGLTHSYPDPYFLACLLFXPENKEVDQDSSLIEKYVSSLNRSFRRQYKHMCRSRQPSTLFXLGQKKGLNSLVHKAEIERYFSEVQDSNSFWQNGVVWEKGEVKFLLCLLDGXAEDKQISLEYGTEAKIKIPVTSVYSGPLRSGRNIERASFYLGFFIEGPLAYGIKVI</TSeq_sequence>
    </TSeq>

    </TSeqSet>
    """
    command_format = 'curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id={}&amp;rettype=fasta&amp;retmode=xml"'
    for prot in proteins:
        if prot in prot2taxid:
            continue
        command = command_format.format(prot)
        cmd = shlex.split(command)
        #try:
        out = subprocess.check_output(cmd)
        #except:
            #print("Something went wrong running command {}".format(command))
            #sys.exit(1)
        res = out.decode()
        taxid, spname = parse_results(res)
        taxid2prot.setdefault(taxid, list()).append(prot)
        taxid2sp[taxid] = spname
        prot2taxid.setdefault(prot, list()).append(taxid)
        time.sleep(0.1)
        
def read_known_prot2taxid(path):
    taxid2prot, taxid2sp, prot2taxid = dict(), dict(), dict()
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            prot = tmp[0]
            taxid = tmp[1]
            spname = tmp[2]
            taxid2prot.setdefault(taxid, list()).append(prot)
            taxid2sp[taxid] = spname
            prot2taxid.setdefault(prot, list()).append(taxid)
    return taxid2prot, taxid2sp, prot2taxid

def replace_tree_id2sp(tree, taxid2sp):
    for taxid in taxid2sp:
        sp = taxid2sp[taxid]
        sp = sp.replace("(", "").replace(")", "").replace(",", "").replace(";", "")
        taxid2sp[taxid] = sp
    
    for taxid in taxid2sp:
    #    #print(taxid)
        node = tree.search_nodes(name=taxid)
        if node != []:
            node[0].name = taxid2sp[taxid]
        else:
            print("Unable to find node for taxid {}".format(taxid))
    return tree

def combine_features(tree, taxid2prot, taxid2sp, prot2taxid, proteins):
    """ create a ete3 tree with protein from jackhmmer
    """ 
    sp2prot = dict()
    for taxid in taxid2sp:
        sp = taxid2sp[taxid]
        sp = sp.replace("(", "").replace(")", "").replace(",", "").replace(";", "")
        taxid2sp[taxid] = sp
        for prot in taxid2prot[taxid]:
            sp2prot.setdefault(sp, list()).append(prot)
            
    # expand taxonomic tree by the number of sequences in each taxa
    for node in tree:
        if node.is_leaf():
            node_sp = node.name
            if node_sp in sp2prot:
                targets = list()
                for prot in sp2prot[node_sp]:
                    for target in proteins[prot]:
                        targets.append(sp+" | "+target)
                subtree = ete3.PhyloTree("({});".format(", ".join(targets)))
                node.add_child(subtree)
    return tree


def main():
    params = get_cmd()
    proteins, fasta = read_inputfile(params.inputfile)
    
    taxid2prot, taxid2sp, prot2taxid = dict(), dict(), dict()
    if params.knownprot2taxid:
        taxid2prot, taxid2sp, prot2taxid = read_known_prot2taxid(params.knownprot2taxid)
        
    retrieve_taxid(proteins, taxid2prot, taxid2sp, prot2taxid)
    
    with open(params.outtaxid, "w") as outf:
        for prot in prot2taxid:
            for taxid in prot2taxid[prot]:
                outf.write("{} {} {}\n".format(prot, taxid, taxid2sp[taxid]))
    
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(list(taxid2prot.keys()), intermediate_nodes=True)


    # check that taxid in tree
    for taxid in taxid2prot:
        node = tree.search_nodes(name=taxid)
        if len(node) < 1:
            print("Unable to find taxonomic id {} in tree".format(taxid), file=sys.stderr)
            sys.exit(1)
    
    # create tree with sp name
    tree = replace_tree_id2sp(tree, taxid2sp)
    tree.render(params.sptree, dpi=300, w=8000)

    # create tree with alignment
    new_tree = combine_features(tree, taxid2prot, taxid2sp, prot2taxid, proteins)
    new_tree.link_to_alignment(fasta, alg_format="fasta")

    
    # clean tree
    to_remove = list()
    for node in tree.traverse():
        if not node.is_root() and not node.is_leaf():
            if len(node.children) == 1:
                to_remove.append(node)
    while to_remove != []:
        for node in to_remove:
            node.delete()
        to_remove = list()
        for node in tree.traverse():
            if not node.is_root() and not node.is_leaf():
                if len(node.children) == 1:
                    to_remove.append(node)

    
    most_distant_leaf, tree_length = new_tree.get_farthest_leaf()
    current_dist = 0
    for postorder, node in new_tree.iter_prepostorder():
        if postorder:
            current_dist -= node.dist
        else:
            if node.is_leaf():
                node.dist += tree_length - (current_dist + node.dist)
            elif node.up: # node is internal
                current_dist += node.dist
    
    ts = TreeStyle()    
    ts.branch_vertical_margin = 10
    ts.allow_face_overlap = False
    ts.show_scale = False
    ts.show_leaf_name = False
    
    new_tree.render(params.alitree, tree_style=ts, dpi=300, w=8000)
    sys.exit(0)

if __name__ == "__main__":
    main()
