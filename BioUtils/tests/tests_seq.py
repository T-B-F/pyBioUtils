from __future__ import division, absolute_import, print_function

import os, sys
import numpy as np
import unittest

from BioUtils import BUseq

class TestBUseq(unittest.TestCase):

    def setup(self):
        self.bad_sequence = "IGEGHPWC-RLCMQN?HCGAD!RHGR*HCQMM."
        self.sequence = "IGEGHPWCRLCMQNHCGADRHGRHCQMM" 
        self.msa_sequences = ["------WCRLCMQNHC----T-HGR-H-LCMQN---",
                              "R----LIMV-L-I-V----RGTY"]
        
    def test_compute_pos_msa2seq(self):
        """ check convertion positions between seq and msa
        """
        self.setup()
        for seq in self.msa_sequences:
            seq_ = seq.replace("-", "")
            for i in range(10):
                start = np.random.randint(0, int(len(seq)/2))
                stop = np.random.randint(int(len(seq)/2)+1, len(seq))
                t_seq = seq_[start: stop]
                new_start, new_stop = BUseq.compute_pos_msa2seq(seq_, start, stop)
                t_seq_ = seq_[new_start: new_stop].replace("-", "")
                self.assertEqual(t_seq_, t_seq)

    def test_compute_pos_seq2msa(self):
        """ check convertion positions between msa and seq
        """
        self.setup()
        for seq in self.msa_sequences:
            seq_ = seq.replace("-", "")
            for i in range(10):
                start = np.random.randint(0, int(len(seq_)/2))
                stop = np.random.randint(int(len(seq_)/2)+1, len(seq_))
                t_seq_ = seq_[start: stop]
                new_start, new_stop = BUseq.compute_pos_seq2msa(seq, start, stop)
                t_seq = seq[new_start: new_stop].replace("-", "")
                self.assertEqual(t_seq_, t_seq)

    def test_compute_offset_pos(self):
        """ check positions
        """
        self.setup()
        seq = self.msa_sequences[0]
        pos = BUseq.compute_offset_pos(seq, 0)
        self.assertEqual(("W", 6), (seq[pos], pos))
        pos = BUseq.compute_offset_pos(seq, 10)
        self.assertEqual(("T", 20), (seq[pos], pos))
        pos = BUseq.compute_offset_pos(seq, 14)
        self.assertEqual(("H", 26), (seq[pos], pos))
        pos = BUseq.compute_offset_pos(seq, 19)
        self.assertEqual(("N", 32), (seq[pos], pos))

        seq = self.msa_sequences[1]
        pos = BUseq.compute_offset_pos(seq, 0)
        self.assertEqual(("R", 0), (seq[pos], pos))
        pos = BUseq.compute_offset_pos(seq, 6)
        self.assertEqual(("I", 12), (seq[pos], pos))
        pos = BUseq.compute_offset_pos(seq, 11)
        self.assertEqual(("Y", 22), (seq[pos], pos))

    def test_transform_seq(self):
        """ check bad characters in MSA are remove """
        self.setup()
        new_seq = BUseq.transform_seq(self.bad_sequence)
        self.assertEqual(new_seq, self.sequence)


if __name__ == "__main__":
    unittest.main()
        
