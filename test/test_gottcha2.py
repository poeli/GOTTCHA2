import unittest
import os
import gottcha.scripts.gottcha2 as gottcha2

class TestGottcha2(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.dirname(os.path.abspath(__file__))
        
    def test_merge_ranges(self):
        ranges = [(1,5), (2,6), (8,10), (11,12)]
        expected = [(1,6), (8,12)]
        result = gottcha2.merge_ranges(ranges)
        self.assertEqual(result, expected)

    def test_parse_sam(self):
        sam_line = "read1\t0\tABC|1|100|12345\t11\t0\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0"
        ref, region, nm, rname, rseq, rq, flag, cigr, pri_flag, valid_flag = gottcha2.parse(sam_line, 0.5)
        self.assertEqual(ref, "ABC|1|100|12345")
        self.assertEqual(region, (11, 20))
        self.assertEqual(nm, 0)
        self.assertTrue(pri_flag)
        self.assertTrue(valid_flag)

    def test_seqReverseComplement(self):
        seq = "ATCG"
        result = gottcha2.seqReverseComplement(seq)
        self.assertEqual(result, "CGAT")

    def test_pile_lvl_zscore(self):
        zscore = gottcha2.pile_lvl_zscore(100, 1000, 100)
        self.assertIsInstance(zscore, float)

if __name__ == '__main__':
    unittest.main()
