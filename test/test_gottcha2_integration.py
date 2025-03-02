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
import src.gottcha2.gottcha2 as gottcha2

class TestGottcha2Integration(unittest.TestCase):
    """Integration tests for GOTTCHA2."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')
        
        # Create a simple mock SAM file
        self.sam_file = os.path.join(self.test_dir, "test.sam")
        with open(self.sam_file, 'w') as f:
            f.write("read1\t0\tABC|1|100|12345\t11\t0\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\n")
            f.write("read2\t0\tABC|1|100|12345\t21\t0\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\n")
            f.write("read3\t0\tXYZ|1|100|67890\t31\t0\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\n")
    
    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.test_dir)
    
    @patch('src.gottcha2.gottcha2.Pool')
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
    
    def test_is_descendant(self):
        """Test checking if a taxid is descendant of another."""
        with patch('src.gottcha2.gottcha2.gt.taxid2fullLineage', return_value='|1|2|1234|5678|'):
            # Test where target is in lineage
            self.assertTrue(gottcha2.is_descendant('5678', '1234'))
            # Test where target is not in lineage
            self.assertFalse(gottcha2.is_descendant('5678', '9999'))

    @patch('src.gottcha2.gottcha2.pd.DataFrame')
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

if __name__ == '__main__':
    unittest.main()
