from __future__ import division, absolute_import, print_function

import os, sys
import numpy as np
import unittest

from BioUtils import BUio

class TestBUseq(unittest.TestCase):

    def setUp(self):
        self.targets = [("fasta/msa_test1.fa", "fasta/seq_test1.fa", 3, 102),
                        ("fasta/msa_test2.fa.gz", "fasta/seq_test2.fa", 3, 102)]

                

    def test_msaio(self):
        """ check msa reading
        """
        # check full reading
        for msafilename, seqfilename, nb_seq, len_seq in self.targets:
            msa_fasta = BUio.read_fastabio(msafilename)
            proteins = list(msa_fasta.keys())
            self.assertEqual(nb_seq, len(msa_fasta))
            self.assertEqual(len(msa_fasta[proteins[0]]), len_seq)

            flat_msa_fasta = BUio.msa2flat(msa_fasta)
            seq_fasta = BUio.read_fastabio(seqfilename)
            seqs_flat = [str(flat_msa_fasta[prot].seq) for prot in proteins]
            seqs = [str(seq_fasta[prot].seq) for prot in proteins]
            self.assertEqual(seqs, seqs_flat)

        # check iteration loading
        for msafilename, seqfilename, nb_seq, len_seq in self.targets:
            seq_fasta = BUio.read_fastabio(msafilename)
            seqs_msa = list()
            seqs_seq = list()
            for record in BUio.read_fastabioit(msafilename):
                seq = str(record.seq)
                seqs_msa.append(seq)
                seqs_seq.append(str(seq_fasta[record.id].seq))
            self.assertEqual(seqs_msa, seqs_seq)



if __name__ == "__main__":
    unittest.main()
        
