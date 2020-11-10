from __future__ import division, absolute_import, print_function

import os, sys
import numpy as np
import unittest

from BioUtils import BUresidue

class TestBUresidue(unittest.TestCase):

    def setUp(self):
        self.AA1 = set(['A', 'C', 'D', 'E', 'F',
                   'G', 'H', 'I', 'K', 'L', 
                   'M', 'N', 'P', 'Q', 'R', 
                   'S', 'T', 'V', 'W', 'Y'])
        self.a2aaa = {'A': "ALA", 'C': "CYS", 'D': "ASP",
                 'E': "GLU", 'F': "PHE", 'G': "GLY",
                 'H': "HIS", 'I': "ILE", 'K': "LYS",
                 'L': "LEU", 'M': "MET", 'N': "ASN",
                 'P': "PRO", 'Q': "GLN", 'R': "ARG",
                 'S': "SER", 'T': "THR", 'V': "VAL",
                 'W': "TRP", 'Y': "TYR"}


    def test_AA(self):
        """ check all common AA are present
        """
        for aa in BUresidue.AA1("classic"):
            self.assertIn(aa, self.AA1)

        for a in BUresidue.AA1("classic"):
            aaa = BUresidue.AA1toAA3(a)
            self.assertEqual(aaa, self.a2aaa[a])

    def test_hydrophobicity(self):
        """ check all aa in hydrophobicity scale and AA1 is equal to AA3
        """
        for scale in ["interface", "octanol", "octanol_interface"]:
            hscale = BUresidue.hydrophobicity_scale(scale)
            for a in self.AA1:
                aaa = self.a2aaa[a]
                self.assertIn(a, hscale)
                self.assertIn(aaa, hscale)
                if a in hscale and aaa in hscale:
                    v = hscale[a]
                    vvv = hscale[a]
                    self.assertEqual(v, vvv)

if __name__ == "__main__":
    unittest.main()
        
