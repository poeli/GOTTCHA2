import unittest
import os
import sys
import tempfile
import shutil
import pandas as pd
from unittest.mock import patch, MagicMock
from pathlib import Path

# Get the path to the project root directory
current_dir = Path(__file__).parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import the module
import gottcha.scripts.gottcha2 as gottcha2

class TestGottcha2Integration(unittest.TestCase):
    """Integration tests for GOTTCHA2."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')
        
        # Create a mock SAM file with proper SAM format
        self.sam_file = os.path.join(self.test_dir, "test.sam")
        with open(self.sam_file, 'w') as f:
            # Header line (required for SAM format)
            f.write("@HD\tVN:1.0\n")
            # Sequence dictionary lines
            f.write("@SQ\tSN:ABC|1|100|12345\tLN:100\n")
            f.write("@SQ\tSN:XYZ|1|100|67890\tLN:100\n")
            # Actual alignment records
            f.write("read1\t0\tABC|1|100|12345\t11\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
            f.write("read2\t0\tABC|1|100|12345\t21\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
            f.write("read3\t0\tXYZ|1|100|67890\t31\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
    
    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.test_dir)
    
    @patch('gottcha.scripts.gottcha2.Pool')
    def test_process_sam_file(self, mock_pool):
        """Test processing of SAM file."""
        # Mock the Pool and its return values
        mock_job = MagicMock()
        mock_pool.return_value.apply_async.return_value = mock_job
        
        # Set up mock data to return from worker
        mock_result = {
            'ABC|1|100|12345': {
                'REGIONS': [(11, 20), (21, 30)],
                'MB': 20,
                'MR': 2,
                'NM': 0
            },
            'XYZ|1|100|67890': {
                'REGIONS': [(31, 40)],
                'MB': 10,
                'MR': 1,
                'NM': 0
            }
        }
        mock_job.get.return_value = mock_result
        
        # Create an actual test file with gottcha2 global variables
        gottcha2.argvs = MagicMock()
        gottcha2.argvs.debug = False
        gottcha2.argvs.silent = True
        gottcha2.begin_t = 0
        gottcha2.logfile = os.path.join(self.test_dir, "test.log")
        
        # Run the function
        result, mapped_reads = gottcha2.process_sam_file(self.sam_file, 1, 0.5)
        
        # Check the results
        self.assertEqual(mapped_reads, 3)  # 2 + 1 from the mocked data
        self.assertEqual(len(result), 2)   # Two reference sequences
        self.assertEqual(result['ABC|1|100|12345']['LL'], 20)  # Linear length of merged regions
        self.assertEqual(result['XYZ|1|100|67890']['LL'], 10)
    
    def test_isin_target_taxa(self):
        """Test checking if a taxid is descendant of another."""
        with patch('gottcha.scripts.gottcha2.gt.taxid2fullLineage', return_value='|1|2|1234|5678|'):
            # Test where target is in lineage
            self.assertTrue(gottcha2.isin_target_taxa('5678', '1234'))
            # Test where target is not in lineage
            self.assertFalse(gottcha2.isin_target_taxa('5678', '9999'))
            # Test where target is in lineage
            self.assertTrue(gottcha2.isin_target_taxa('5678', ['1234','9999']))


    @patch('gottcha.scripts.gottcha2.pd.DataFrame')
    def test_group_refs_to_strains(self, mock_df):
        """Test grouping references by strain."""
        # Mock the pandas operations
        mock_df.from_dict.return_value = mock_df
        mock_df.rename.return_value = mock_df
        mock_df.str.rstrip.return_value = mock_df
        mock_df.str.split.return_value.expand = True
        
        # Set up a mock database stats
        gottcha2.db_stats = {'12345': 1000, '67890': 500}
        
        # Test the function with mock data
        test_data = {
            'ABC|1|100|12345': {'MB': 100, 'MR': 10, 'NM': 2, 'LL': 50},
            'XYZ|1|100|67890': {'MB': 50, 'MR': 5, 'NM': 1, 'LL': 25}
        }
        
        # Call the function (will use mocked pandas operations)
        gottcha2.group_refs_to_strains(test_data)
        
        # Verify the expected pandas calls were made (basic verification)
        mock_df.from_dict.assert_called_once()

    def test_multiple_taxid_extraction(self):
        """Test extracting reads matching multiple taxids."""
        # Mock the taxonomy lineage lookup
        with patch('gottcha.scripts.gottcha2.gt.taxid2fullLineage') as mock_lineage, \
             patch('gottcha.scripts.gottcha2.parse') as mock_parse:
            
            def lineage_side_effect(taxid, space2underscore=False):
                lineages = {
                    '12345': '|1|2|12345|',
                    '67890': '|1|2|67890|',
                    '99999': '|1|2|99999|'
                }
                return lineages.get(taxid, '')
            mock_lineage.side_effect = lineage_side_effect
            
            def parse_side_effect(line, matchFactor):
                if line.startswith('@'):
                    return None, None, None, None, None, None, None, None, None, None
                parts = line.split('\t')
                ref = parts[2]
                region = (int(parts[3]), int(parts[3])+10)  # Assuming 10bp matches
                nm = 0
                rname = parts[0]
                rseq = parts[9]
                rq = parts[10]
                flag = int(parts[1])
                cigr = parts[5]
                return ref, region, nm, rname, rseq, rq, flag, cigr, True, True
            mock_parse.side_effect = parse_side_effect
            
            # Test extraction with multiple taxids using context manager
            with open(self.sam_file, 'r') as f:
                result, read_count = gottcha2.ReadExtractWorker(
                    f.name,
                    0,
                    os.path.getsize(self.sam_file),
                    "12345,67890",
                    0.5
                )
                
                # Parse results
                result_lines = [line for line in result.strip().split('\n') if line]
                reads = [line for line in result_lines if line.startswith('@')]
                
                # Verify correct reads were extracted
                self.assertEqual(read_count, 3, "Should extract 3 reads")
                
                # Check taxids in extracted reads
                extracted_taxids = set()
                for read in reads:
                    if '12345' in read:
                        extracted_taxids.add('12345')
                    elif '67890' in read:
                        extracted_taxids.add('67890')
                        
                self.assertEqual(
                    extracted_taxids,
                    {'12345', '67890'},
                    "Should extract reads from both specified taxids"
                )
                
                # Verify exclusion of unspecified taxid
                self.assertTrue(
                    all('99999' not in line for line in result_lines),
                    "Should not extract reads from unspecified taxid 99999"
                )

    def test_taxid_file_input(self):
        """Test extracting reads using taxids from a file."""
        # Create a temporary taxid file
        taxid_file = os.path.join(self.test_dir, "taxids.txt")
        with open(taxid_file, 'w') as f:
            f.write("12345\n")
            f.write("67890\n")
            f.write("\n")  # Empty line to test handling
            f.write("# Comment to ignore\n")
        
        # Test extraction using file input
        taxids = gottcha2.parse_taxids(f'@{taxid_file}')
                        
        self.assertEqual(len(taxids), 2, f"Should have 2 taxids")
        self.assertEqual(
            taxids,
            ['12345', '67890'],
            "Should extract taxids listed in file"
        )

if __name__ == '__main__':
    unittest.main()
