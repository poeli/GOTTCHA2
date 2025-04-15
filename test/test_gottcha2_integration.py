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
        mock_result = ({
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
        }, 3, 0, 0)  # Return tuple matching updated function signature
        
        mock_job.get.return_value = mock_result
        
        # Create an actual test file with gottcha2 global variables
        gottcha2.argvs = MagicMock()
        gottcha2.argvs.debug = False
        gottcha2.argvs.silent = True
        gottcha2.begin_t = 0
        gottcha2.logfile = os.path.join(self.test_dir, "test.log")
        
        # Run the function
        result, mapped_reads, tol_alignment_count, tol_invalid_match_count, tol_exclude_acc_count = gottcha2.process_sam_file(self.sam_file, 1, 0.5)
        
        # Check the results
        self.assertEqual(mapped_reads, 3)  # 2 + 1 from the mocked data
        self.assertEqual(len(result), 2)   # Two reference sequences
        self.assertEqual(result['ABC|1|100|12345']['LL'], 20)  # Linear length of merged regions
        self.assertEqual(result['XYZ|1|100|67890']['LL'], 10)
        self.assertEqual(tol_alignment_count, 3)
        self.assertEqual(tol_invalid_match_count, 0)
        self.assertEqual(tol_exclude_acc_count, 0)
    
    def test_isin_target_taxa(self):
        """Test checking if a taxid is descendant of another."""
        with patch('gottcha.scripts.gottcha2.gt.taxid2fullLineage', return_value='|1|2|1234|5678|'):
            # Test where target is in lineage
            self.assertTrue(gottcha2.isin_target_taxa('5678', ['1234']))
            # Test where target is not in lineage
            self.assertFalse(gottcha2.isin_target_taxa('5678', ['9999']))
            # Test where target is in lineage
            self.assertTrue(gottcha2.isin_target_taxa('5678', ['1234','9999']))

    @patch('gottcha.scripts.gottcha2.pd.DataFrame')
    @patch('gottcha.scripts.gottcha2.logging.fatal')
    @patch('gottcha.scripts.gottcha2.sys.exit')
    def test_group_refs_to_strains(self, mock_exit, mock_fatal, mock_df):
        """Test grouping references by strain."""
        # Mock the pandas operations
        mock_df.from_dict.return_value = mock_df
        mock_df.rename.return_value = mock_df
        mock_df.str.rstrip.return_value = mock_df
        mock_df.str.split.return_value.expand = True
        
        # Create a mock Series for eq(0).any() check to return False
        mock_series = MagicMock()
        mock_series.eq.return_value.any.return_value = False
        mock_df.__getitem__.return_value = mock_series
        
        # Set up a mock database stats with the indexes set properly
        gottcha2.df_stats = pd.DataFrame([
            ['species', '12345', 1000, 10000],
            ['species', '67890', 500, 5000]
        ], columns=['DB_level', 'Taxid', 'TotalLength', 'GenomeSize'])
        gottcha2.df_stats.set_index('Taxid', inplace=True)
        
        # Test the function with mock data
        test_data = {
            'ABC|1|100|12345': {'MB': 100, 'MR': 10, 'NM': 2, 'LL': 50},
            'XYZ|1|100|67890': {'MB': 50, 'MR': 5, 'NM': 1, 'LL': 25}
        }
        
        # Call the function (will use mocked pandas operations)
        gottcha2.group_refs_to_strains(test_data)
        
        # Verify the expected pandas calls were made
        mock_df.from_dict.assert_called_once()
        # mock_series.eq.assert_called_once_with(0)
        
        # # Verify sys.exit was not called
        # mock_exit.assert_not_called()
        # mock_fatal.assert_not_called()

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
            
            def parse_side_effect(line, matchFactor, excluded_acc_list=None):
                if line.startswith('@'):
                    return None, None, None, None, None, None, None, None, None, None, None
                parts = line.split('\t')
                ref = parts[2]
                region = (int(parts[3]), int(parts[3])+10)  # Assuming 10bp matches
                nm = 0
                rname = parts[0]
                rseq = parts[9]
                rq = parts[10]
                flag = int(parts[1])
                cigr = parts[5]
                return ref, region, nm, rname, rseq, rq, flag, cigr, True, True, True
            mock_parse.side_effect = parse_side_effect
            
            # Test extraction with multiple taxids using context manager
            with open(self.sam_file, 'r') as f:
                result, read_count = gottcha2.ReadExtractWorker(
                    f.name,
                    0,
                    os.path.getsize(self.sam_file),
                    ["12345", "67890"],  # Pass as a list
                    0.5,
                    None
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
        
        # Create mock taxonomy file data
        mock_taxa_df = pd.DataFrame({
            'LEVEL': ['species', 'species'],
            'NAME': ['Species1', 'Species2'],
            'TAXID': ['12345', '67890'],
            'NOTE': ['', '']
        })

        # Mock logging and file operations
        with patch('gottcha.scripts.gottcha2.print_message') as mock_print, \
             patch('gottcha.scripts.gottcha2.pd.read_csv', return_value=mock_taxa_df):
            
            # Test with empty DataFrame and file input
            gottcha2.argvs = MagicMock()
            gottcha2.argvs.silent = True
            gottcha2.begin_t = 0
            gottcha2.logfile = os.path.join(self.test_dir, "test.log")
            
            empty_df = pd.DataFrame()
            taxa_dict, qualified_taxids = gottcha2.parse_taxids(f'@{taxid_file}', empty_df, "mock_full.tsv")
            
            # Verify the results
            self.assertEqual(len(qualified_taxids), 2)
            self.assertListEqual(sorted(qualified_taxids), ['12345', '67890'])
            self.assertEqual(len(taxa_dict), 2)
            
            # Verify taxonomy dictionary contents
            self.assertEqual(taxa_dict['12345']['level'], 'species')
            self.assertEqual(taxa_dict['12345']['name'], 'Species1')
            self.assertEqual(taxa_dict['67890']['level'], 'species')
            self.assertEqual(taxa_dict['67890']['name'], 'Species2')

            # Test with pre-loaded DataFrame
            taxa_dict2, qualified_taxids2 = gottcha2.parse_taxids(f'@{taxid_file}', mock_taxa_df, "mock_full.tsv")
            self.assertEqual(len(qualified_taxids2), 2)
            self.assertListEqual(sorted(qualified_taxids2), ['12345', '67890'])
            
            # Test with direct taxid input
            taxa_dict3, qualified_taxids3 = gottcha2.parse_taxids('12345,67890', mock_taxa_df, "mock_full.tsv")
            self.assertEqual(len(qualified_taxids3), 2)
            self.assertListEqual(sorted(qualified_taxids3), ['12345', '67890'])
            
            # Test with 'all' option
            taxa_dict4, qualified_taxids4 = gottcha2.parse_taxids('all', mock_taxa_df, "mock_full.tsv")
            self.assertEqual(len(qualified_taxids4), 2)
            self.assertListEqual(sorted(qualified_taxids4), ['12345', '67890'])

    def test_load_excluded_acc_list(self):
        """Test loading accession list to exclude."""
        # Create a temporary file with accessions
        acc_file = os.path.join(self.test_dir, "exclude_acc.txt")
        with open(acc_file, 'w') as f:
            f.write("ABC\nXYZ\nPQR\n")
        
        # Test loading the file
        with open(acc_file) as f:
            exclude_set = gottcha2.load_excluded_acc_list(f)
        
        self.assertEqual(len(exclude_set), 3)
        self.assertIn("ABC", exclude_set)
        self.assertIn("XYZ", exclude_set)
        self.assertIn("PQR", exclude_set)
        
        # Test loading empty file
        empty_file = os.path.join(self.test_dir, "empty.txt")
        with open(empty_file, 'w') as f:
            pass
            
        with open(empty_file) as f:
            with self.assertLogs(level='WARNING'):
                empty_set = gottcha2.load_excluded_acc_list(f)
                
        self.assertEqual(len(empty_set), 0)

    def test_exclude_accessions_in_parse(self):
        """Test that parse() correctly excludes accessions."""
        # Create SAM line with accession ABC
        sam_line = "read1\t0\tABC|1|100|12345\t11\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0"
        
        # Test without exclusion
        _, _, _, _, _, _, _, _, _, valid_match_flag1, valid_acc_flag1 = gottcha2.parse(sam_line, 0.5)
        self.assertTrue(valid_match_flag1)
        self.assertTrue(valid_acc_flag1)
        
        # Test with exclusion
        exclude_set = {"ABC"}
        _, _, _, _, _, _, _, _, _, valid_match_flag2, valid_acc_flag2 = gottcha2.parse(sam_line, 0.5, exclude_set)
        self.assertTrue(valid_match_flag2)
        self.assertFalse(valid_acc_flag2)
        
        # Test with a different accession that's not excluded
        sam_line2 = "read3\t0\tDEF|1|100|67890\t11\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0"
        _, _, _, _, _, _, _, _, _, valid_match_flag3, valid_acc_flag3 = gottcha2.parse(sam_line2, 0.5, exclude_set)
        self.assertTrue(valid_match_flag3)
        self.assertTrue(valid_acc_flag3)

    def test_worker_with_excluded_acc_list(self):
        """Test that worker() correctly excludes accessions."""
        # Create a temporary SAM file with multiple accessions
        test_sam = os.path.join(self.test_dir, "acc_test.sam")
        with open(test_sam, 'w') as f:
            f.write("read1\t0\tABC|1|100|12345\t11\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
            f.write("read2\t0\tXYZ|1|100|67890\t21\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
            f.write("read3\t0\tDEF|1|100|54321\t31\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\tNM:i:0\tAS:i:10\tXS:i:0\n")
        
        # Test without exclusion
        result_no_exclude, lines_count1, invalid_match_count1, exclude_acc_count1 = gottcha2.worker(test_sam, 0, os.path.getsize(test_sam), 0.5)
        self.assertEqual(len(result_no_exclude), 3)
        self.assertIn("ABC|1|100|12345", result_no_exclude)
        self.assertIn("XYZ|1|100|67890", result_no_exclude)
        self.assertIn("DEF|1|100|54321", result_no_exclude)
        self.assertEqual(lines_count1, 3)
        self.assertEqual(invalid_match_count1, 0)
        self.assertEqual(exclude_acc_count1, 0)
        
        # Test with exclusion
        exclude_set = {"ABC", "XYZ"}
        result_with_exclude, lines_count2, invalid_match_count2, exclude_acc_count2 = gottcha2.worker(test_sam, 0, os.path.getsize(test_sam), 0.5, exclude_set)
        self.assertEqual(len(result_with_exclude), 1)
        self.assertNotIn("ABC|1|100|12345", result_with_exclude)
        self.assertNotIn("XYZ|1|100|67890", result_with_exclude)
        self.assertIn("DEF|1|100|54321", result_with_exclude)
        self.assertEqual(lines_count2, 3)
        self.assertEqual(invalid_match_count2, 0)
        self.assertEqual(exclude_acc_count2, 2)  # Two accessions excluded

    @patch('gottcha.scripts.gottcha2.Pool')
    @patch('gottcha.scripts.gottcha2.print_message')
    def test_process_sam_with_exclude_list(self, mock_print, mock_pool):
        """Test processing SAM file with exclude accession list."""
        # Set up mocks
        mock_job = MagicMock()
        mock_pool.return_value.apply_async.return_value = mock_job
        
        # Create result with one accession excluded - matching the updated return signature
        mock_result = ({
            'DEF|1|100|54321': {
                'REGIONS': [(31, 40)],
                'MB': 10,
                'MR': 1,
                'NM': 0
            }
        }, 3, 0, 2)  # Third value is now exclude_acc_count
        
        mock_job.get.return_value = mock_result
        
        # Create a temporary accession list file
        acc_file = os.path.join(self.test_dir, "exclude_acc.txt")
        with open(acc_file, 'w') as f:
            f.write("ABC\nXYZ\n")
        
        # Set up gottcha2 global variables
        gottcha2.argvs = MagicMock()
        gottcha2.argvs.debug = False
        gottcha2.argvs.silent = True
        gottcha2.begin_t = 0
        gottcha2.logfile = os.path.join(self.test_dir, "test.log")
        
        # Load exclusion list and run the function
        with open(acc_file) as f:
            exclude_set = gottcha2.load_excluded_acc_list(f)
            
        result, mapped_reads, tol_alignment_count, tol_invalid_match_count, tol_exclude_acc_count = gottcha2.process_sam_file(
            self.sam_file, 1, 0.5, exclude_set
        )
        
        # Check the results
        mock_pool.return_value.apply_async.assert_called_with(
            gottcha2.worker, 
            (self.sam_file, mock_pool.return_value.apply_async.call_args[0][1][1], mock_pool.return_value.apply_async.call_args[0][1][2], 0.5, exclude_set)
        )
        
        # Only one reference should be in the result
        self.assertEqual(len(result), 1)
        self.assertIn('DEF|1|100|54321', result)
        self.assertEqual(tol_alignment_count, 3)
        self.assertEqual(tol_invalid_match_count, 0)
        self.assertEqual(tol_exclude_acc_count, 2)  # Two accessions excluded

if __name__ == '__main__':
    unittest.main()
