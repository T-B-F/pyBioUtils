#!/usr/bin/env python

import sys
import numpy as np

def AA1toAA3(aa, verbose=True):
    if(aa == 'A'):
        return "ALA"
    elif(aa == 'C'):
        return "CYS"
    elif(aa == 'D'):
        return "ASP"
    elif(aa == 'E'):
        return("GLU")
    elif(aa == 'F') :
        return("PHE")
    elif(aa == 'G') :
        return("GLY")
    elif(aa == 'H') :
        return("HIS")
    elif(aa == 'I') :
        return("ILE")
    elif(aa == 'K') :
        return("LYS")
    elif(aa == 'L') :
        return("LEU")
    elif(aa == 'M') :
        return("MET")
    elif(aa == 'N') :
        return("ASN")
    elif(aa == 'P') :
        return("PRO")
    elif(aa == 'Q') :
        return("GLN")
    elif(aa == 'R') :
        return("ARG")
    elif(aa == 'S') :
        return("SER")
    elif(aa == 'T') :
        return("THR")
    elif(aa == 'V') :
        return("VAL")
    elif(aa == 'W') :
        return("TRP")
    elif(aa == 'Y') :
        return("TYR")
    else :
        if verbose :
            print("The amino acid {} is unknown!".format(aa), file=sys.stderr)
        return "UNK"

def AA3toAA1(aa, verbose=True):
    if(aa == 'ALA'):
        return "A"
    elif(aa == 'CYS'):
        return "C"
    elif(aa == 'ASP'):
        return "D"
    elif(aa == 'GLU'):
        return("E")
    elif(aa == 'PHE') :
        return("F")
    elif(aa == 'GLY') :
        return("G")
    elif(aa == 'HIS') :
        return("H")
    elif(aa == 'ILE') :
        return("I")
    elif(aa == 'LYS') :
        return("K")
    elif(aa == 'LEU') :
        return("L")
    elif(aa == 'MET') :
        return("M")
    elif(aa == 'ASN') :
        return("N")
    elif(aa == 'PRO') :
        return("P")
    elif(aa == 'GLN') :
        return("Q")
    elif(aa == 'ARG') :
        return("R")
    elif(aa == 'SER') :
        return("S")
    elif(aa == 'THR') :
        return("T")
    elif(aa == 'VAL') :
        return("V")
    elif(aa == 'TRP') :
        return("W")
    elif(aa == 'TYR') :
        return("Y")
    # non standart amino acid  / translation accorded to ProDyn web site http://www.csb.pitt.edu/prody/reference/atomic/flags.html
    elif (aa == 'ASX' ) : # asparagine or aspartic acid
        return "B"     
    elif (aa == 'GLX' ) : # glutamine or glutamic acid
        return "Z"     
    elif (aa == 'CSO') :  # S-hydroxycysteine
        return "C"     
    elif (aa == 'HIP' ) : # ND1-phosphohistidine
        return "H"     
    elif (aa == 'HSD' ) : # prototropic tautomer of histidine, H on ND1 (CHARMM)
        return "H"     
    elif (aa == 'HSE' ) : # prototropic tautomer of histidine, H on NE2 (CHARMM)
        return "H"     
    elif (aa == 'HSP' ) : # protonated histidine
        return "H"     
    elif (aa == 'MSE' ) : # selenomethionine
        return "X" 
    elif (aa == 'SEC' ) : # selenocysteine
        return "X"     
    elif (aa == 'SEP' ) : # phosphoserine
        return "S"     
    elif (aa == 'TPO' ) : # phosphothreonine
        return "T"     
    elif (aa == 'PTR' ) : # O-phosphotyrosine
        return "Y"     
    elif (aa == 'XLE' ) : #  leucine or isoleucine
        return "L"    # or I
    elif (aa == 'XAA') : #  unspecified or unknown
        return "X"
    elif (aa == 'UNK'):
        return "X"
    else :
        if verbose :
            print("The amino acid {} is unknown!".format(aa), file=sys.stderr)
        return "X"
        #raise ValueError( "The amino acid ",aa," is unknown!" )        


veryHydrophobic = ["VAL","ILE", "LEU", "MET", "PHE", "TRP", "CYS"]
lessHydrophobic = ["ALA", "TYR", "HIS", "THR", "SER", "PRO", "GLY"]
hydrophobic = veryHydrophobic+lessHydrophobic
polar = ["ARG", "LYS", "ASP", "GLU", "ASN", "GLN"]
lessPolar = ["HIS", "ALA", "TYR", "THR", "SER", "PRO", "GLY" ]
tiny = ["GLY", "ALA", "SER", "PRO"]
small = ["THR","ASP","ASN"]
aliphatic = ["ILE", "VAL", "LEU", "ALA", "PRO" ]
aliphaticExtended = aliphatic + ["MET"]
chargedNeg = [ "GLU", "ASP"]
chargedPlus = ["ARG", "LYS"]
hydrophilic = ['GLN', 'ASN', 'LYS', 'ASP', 'ARG', 'GLU']
chargedPlusExtended = chargedPlus+["HIS"]
aromatic = ["PHE", "TRP", "TYR", "HIS" ]

_AA3 = {"all": ['CYS', 'GLN', 'ILE', 'SER', 'VAL', 'GLY', 'ASN', 'PRO', 'LYS', 'ASP', 'THR', 'PHE', 'ALA', 'MET', 'HIS', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR', 'UNK'],
        "classic": ['CYS', 'GLN', 'ILE', 'SER', 'VAL', 'GLY', 'ASN', 'PRO', 'LYS', 'ASP', 'THR', 'PHE', 'ALA', 'MET', 'HIS', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR'],
        }

_AA1 = {"all": [AA3toAA1(aa) for aa in _AA3["all"]],
        "classic": [AA3toAA1(aa) for aa in _AA3["classic"]],
        }

def AA3(mode="all"):
    if mode in _AA3:
        return _AA3[mode]
    else:
        print("Error, mode={} not found. Choose between: {}".format(mode, "/".join(list(_AA3.keys()))), file=sys.stderr)


def AA1(mode="all"):
    if mode in _AA1:
        return _AA1[mode]
    else:
        print("Error, mode={} not found. Choose between: {}".format(mode, "/".join(list(_AA1.keys()))), file=sys.stderr)

dicAA1Type = {
  "hydrophobic" : ["V","I", "L", "M", "F", "W", "C", "Y", "H", "T" ],
  "hydrophilic" : ['Q', 'N', 'K', 'D', 'R', 'E'],
  "tiny" : ["G", "A", "S", "P"] ,
}

dicAA3Type = {
  "hydrophobic" : ["VAL","ILE", "LEU", "MET", "PHE", "TRP", "CYS", "ALA", "TYR", "HIS", "THR", "SER", "PRO", "GLY" ],
  "hydrophilic" : ['GLN', 'ASN', 'LYS', 'ASP', 'ARG', 'GLU'],
  "tiny" : ["GLY", "ALA", "SER", "PRO"] ,
}

dpdb_freq = {"A": 0.0859929227395,
            "C": 0.0206946259759,
            "E": 0.0648133592624,
            "D": 0.0549383608817,
            "G": 0.0830137574641,
            "F": 0.0378506807406,
            "I": 0.0547752297123,
            "H": 0.0260798136952,
            "K": 0.057416225372,
            "M": 0.0230113512204,
            "L": 0.0871762434274,
            "N": 0.0408935875728,
            "Q": 0.0364996548365,
            "P": 0.045826092199,
            "S": 0.0602061283907,
            "R": 0.0504294572634,
            "T": 0.05511062539,
            "W": 0.0131711008412,
            "V": 0.068839479608,
            "Y": 0.0332613034068}

class Frequency(object):
    
    def __init__(self, name=""):
        self._name = name
        if self._name == "uniform":
            self._aa = AA1(mode="classic")
            self._freq = np.array([1/len(self._aa) for aa in self._aa])
            self._freq /= self._freq.sum()
        elif self._name == "pdb":
            self._aa = AA1(mode="classic")
            self._freq = np.array([dpdb_freq[aa] for aa in self._aa])
            self._freq /= self._freq.sum()
        else:
            self._freq = list()
            self._aa = list()
        
    def from_file(self, path):
        self._freq = list()
        self._aa = list()
        with open(path) as inf:
            for line in inf:
                tmp = line.split()
                self._aa.append(tmp[0])
                self._freq.append(float(tmp[1]))
        self._freq = np.array(self._freq)
        self._normalize()
        
    def _normalize(self):
        self._freq /= self._freq.sum()
    
    def todict(self):
        freq = dict()
        for i, aa in enumerate(self._aa):
            freq[aa] = self._freq[i]
        return freq
    
    def from_dict(self, data):
        """ upload frequency from dictionary
        """
        for aa in data:
            self._aa.append(aa)
            self._freq.append(data[aa])
        self._freq = np.array(self._freq)
        self._normalize()
        
    def sequence(self, size, N=1):
        """ generate N random sequence of size "size" using the set of amino acids in aa
        and the corresponding freq 
        """
        if N == 1:
            seq = "".join(np.random.choice(self._aa, p=self._freq) for n in range(size))
        else:
            seq = ["".join(np.random.choice(self._aa, p=self._freq) for n in range(size)) for n in range(N)]
        return seq
    
    def ite_sequence(self, size, N=1):
        """ generate N random sequence of size "size" using the set of amino acids in aa
        and the corresponding freq 
        """
        for i in range(N):
            seq = "".join(np.random.choice(self._aa, p=self._freq) for n in range(size))
            yield seq
            
