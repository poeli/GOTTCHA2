import unittest
import os
import sys
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from pathlib import Path
from io import StringIO
import pandas as pd

# Get the path to the project root directory
current_dir = Path(__file__).parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import the module
import gottcha.scripts.gottcha2 as gottcha2

class TestGottcha2CLI(unittest.TestCase):
    """Test command-line interface functionality of GOTTCHA2."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        # Mock database file
        self.db_path = os.path.join(self.test_dir, "test_db.species.fna")
        # Create a dummy .mmi file to pass database existence check
        with open(f"{self.db_path}.mmi", 'w') as f:
            f.write("dummy minimap2 index")
        # Create a dummy stats file
        with open(f"{self.db_path}.stats", 'w') as f:
            f.write("species\tE. coli\t12345\tBacteria\t1\t100\t50\t1000\t5000000\n")
            f.write("strain\tE. coli str1\t45678\tBacteria\t1\t100\t50\t1000\t5000000\n")
        
        # Create mock input files
        self.test_fastq = os.path.join(self.test_dir, "test.fastq")
        with open(self.test_fastq, 'w') as f:
            f.write("@seq1\nACGT\n+\nHHHH\n")
            
        self.test_sam = os.path.join(self.test_dir, "test.gottcha_species.sam")
        with open(self.test_sam, 'w') as f:
            f.write("@HD\tVN:1.0\n")
            f.write("@SQ\tSN:ABC|1|100|12345\tLN:100\n")
            f.write("read1\t0\tABC|1|100|12345\t11\t60\t5S10M3S\t*\t0\t0\tGGGGGCCCCCCCCCGGG\tHHHHHHHHHHHHHHHHH\n")
            
        # Create mock accession exclusion list
        self.exclusion_list = os.path.join(self.test_dir, "exclude.txt")
        with open(self.exclusion_list, 'w') as f:
            f.write("ABC\nXYZ\n")
    
    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.test_dir)
    
    def test_parse_params(self):
        """Test parameter parsing."""
        # Test version flag
        with patch('sys.exit') as mock_exit:
            gottcha2.parse_params("1.0.0", ['-v'])
            mock_exit.assert_called()
        
        # Test database requirement
        with self.assertRaises(SystemExit):
            gottcha2.parse_params("1.0.0", ['-i', 'test.fastq'])
        
        # Test input and SAM incompatibility
        with self.assertRaises(SystemExit):
            gottcha2.parse_params("1.0.0", ['-i', 'test.fastq', '-s', 'test.sam', '-d', self.db_path])
        
        # Test valid params with input file
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species'])
        self.assertEqual(args.database, self.db_path)
        self.assertEqual(args.dbLevel, 'species')
        
        # Test noCutoff option
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', '-nc'])
        self.assertEqual(args.minCov, 0)
        self.assertEqual(args.minReads, 0)
        self.assertEqual(args.minLen, 0)
        self.assertEqual(args.matchFactor, 0)
        self.assertEqual(args.maxZscore, 0)
        
        # Test nanopore option
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', '-np'])
        self.assertEqual(args.presetx, 'map-ont')
        self.assertEqual(args.minReads, 1)
        self.assertEqual(args.matchFactor, 0)
        
        # Test auto-detection of dbLevel from database name
        db_with_level = os.path.join(self.test_dir, "test_db.species.fna")
        with open(f"{db_with_level.replace('.fna', '')}.mmi", 'w') as f:
            f.write("dummy")
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', db_with_level])
        self.assertEqual(args.dbLevel, 'species')
        
        # Test auto-detection of dbLevel from SAM file name
        args = gottcha2.parse_params("1.0.0", ['-s', self.test_sam, '-d', self.db_path])
        self.assertEqual(args.dbLevel, 'species')
        
        # Test accession exclusion list
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', '-A', self.exclusion_list])
        self.assertEqual(args.accExclusionList.name, self.exclusion_list)
        
        # Test extract option with file
        taxid_file = os.path.join(self.test_dir, "taxids.txt")
        with open(taxid_file, 'w') as f:
            f.write("12345\n")
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', 
                                        '-e', f'@{taxid_file}:10:fasta'])
        self.assertEqual(args.extract, f'@{taxid_file}:10:fasta')
        
        # Test extractFullRef option
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', 
                                        '-ef'])
        self.assertEqual(args.extract, 'all:20:fasta')
        
        # Test extract and extractFullRef incompatibility
        with self.assertRaises(SystemExit):
            gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', 
                                    '-e', '12345', '-ef'])
        
        # Test removeMultipleHits auto mode
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species'])
        self.assertEqual(args.removeMultipleHits, 'yes')
        
        args = gottcha2.parse_params("1.0.0", ['-s', self.test_sam, '-d', self.db_path, '-l', 'species'])
        self.assertEqual(args.removeMultipleHits, 'no')
        
        # Test explicit removeMultipleHits
        args = gottcha2.parse_params("1.0.0", ['-i', self.test_fastq, '-d', self.db_path, '-l', 'species', 
                                        '-rm', 'no'])
        self.assertEqual(args.removeMultipleHits, 'no')
    
    @patch('subprocess.check_call')
    def test_dependency_check(self, mock_check_call):
        """Test dependency checking."""
        # Test successful dependency check
        gottcha2.dependency_check("minimap2")
        mock_check_call.assert_called_once()
        
        # Test failed dependency check
        mock_check_call.side_effect = Exception("Command not found")
        with patch('sys.stderr', new=StringIO()) as mock_stderr:
            with self.assertRaises(SystemExit):
                gottcha2.dependency_check("nonexistent_command")
            self.assertIn("[ERROR]", mock_stderr.getvalue())
    
    @patch('gottcha.scripts.gottcha2.subprocess.Popen')
    def test_readMapping(self, mock_popen):
        """Test read mapping function."""
        # Mock the subprocess.Popen
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("", "")
        mock_process.poll.return_value = 0  # Success exit code
        mock_popen.return_value = mock_process
        
        # Create a file-like object for input reads
        class MockFile:
            def __init__(self, name):
                self.name = name
        
        reads = [MockFile(self.test_fastq)]
        logfile = os.path.join(self.test_dir, "test.log")
        samfile = os.path.join(self.test_dir, "test.sam")
        
        # Test the function
        exit_code, cmd, errs = gottcha2.readMapping(
            reads, self.db_path, 1, 10, 'sr', samfile, logfile
        )
        
        # Check the results
        self.assertEqual(exit_code, 0)
        self.assertIn("minimap2", cmd)
        self.assertIn(self.db_path, cmd)
        mock_popen.assert_called_once()
    
    def test_load_excluded_acc_list(self):
        """Test loading excluded accession list."""
        with open(self.exclusion_list) as f:
            excluded_acc = gottcha2.load_excluded_acc_list(f)
            
        self.assertEqual(len(excluded_acc), 2)
        self.assertIn("ABC", excluded_acc)
        self.assertIn("XYZ", excluded_acc)
        
        # Test with empty file
        empty_file = os.path.join(self.test_dir, "empty.txt")
        with open(empty_file, 'w') as f:
            pass
        
        with open(empty_file) as f:
            with self.assertLogs(level='WARNING'):
                empty_set = gottcha2.load_excluded_acc_list(f)
                
        self.assertEqual(len(empty_set), 0)
    
    def test_loadDatabaseStats(self):
        """Test loading database statistics."""
        # Test with valid stats file
        stats_file = os.path.join(self.test_dir, "valid_stats.tsv")
        with open(stats_file, 'w') as f:
            f.write("species\tE.coli\t12345\tBacteria\t1\t100\t50\t1000\t5000000\t\n")
            f.write("species\tS.aureus\t67890\tBacteria\t1\t100\t50\t2000\t3000000\t\n")

            
        df_stats = gottcha2.loadDatabaseStats(stats_file)
        self.assertEqual(len(df_stats), 2)
        self.assertEqual(df_stats.index[0], '12345')
        self.assertEqual(df_stats.loc['12345', 'TotalLength'], 1000)
        self.assertEqual(df_stats.loc['67890', 'GenomeSize'], 3000000)
        
        # Test with missing genome size column
        stats_file_no_genome = os.path.join(self.test_dir, "no_genome_stats.tsv")
        with open(stats_file_no_genome, 'w') as f:
            f.write("species\tE.coli\t12345\tBacteria\t1\t100\t50\t1000\t\n")
            
        df_stats_no_genome = gottcha2.loadDatabaseStats(stats_file_no_genome)
        self.assertEqual(df_stats_no_genome.loc['12345', 'GenomeSize'], 0)
        
        # Test with non-numeric genome size
        stats_file_bad_genome = os.path.join(self.test_dir, "bad_genome_stats.tsv")
        with open(stats_file_bad_genome, 'w') as f:
            f.write("species\tE.coli\t12345\tBacteria\t1\t100\t50\t1000\tN/A\n")
            
        df_stats_bad_genome = gottcha2.loadDatabaseStats(stats_file_bad_genome)
        self.assertEqual(df_stats_bad_genome.loc['12345', 'GenomeSize'], 0)
    
    def test_parse_taxids(self):
        """Test parsing taxids from command line or file."""
        # Setup mock taxonomy data
        mock_taxa_df = pd.DataFrame({
            'LEVEL': ['species', 'species'],
            'NAME': ['E.coli', 'S.aureus'],
            'TAXID': ['12345', '67890'],
            'NOTE': ['', '']
        })
        
        # Create a taxid file
        taxid_file = os.path.join(self.test_dir, "taxids.txt")
        with open(taxid_file, 'w') as f:
            f.write("12345\n")
            f.write("# Comment line\n")
            f.write("\n")  # Empty line
        
        with patch('gottcha.scripts.gottcha2.print_message') as mock_print:
            gottcha2.argvs = MagicMock()
            gottcha2.argvs.silent = True
            gottcha2.begin_t = 0
            gottcha2.logfile = os.path.join(self.test_dir, "test.log")
            
            # Test file input
            taxa_dict, qualified_taxids = gottcha2.parse_taxids(f'@{taxid_file}', mock_taxa_df, "mock_file.tsv")
            
            self.assertEqual(len(qualified_taxids), 1)
            self.assertEqual(qualified_taxids[0], '12345')
            self.assertEqual(taxa_dict['12345']['level'], 'species')
            self.assertEqual(taxa_dict['12345']['name'], 'E.coli')
            
            # Test comma-separated input
            taxa_dict2, qualified_taxids2 = gottcha2.parse_taxids('12345,67890', mock_taxa_df, "mock_file.tsv")
            
            self.assertEqual(len(qualified_taxids2), 2)
            self.assertIn('12345', qualified_taxids2)
            self.assertIn('67890', qualified_taxids2)
            
            # Test 'all' option
            taxa_dict3, qualified_taxids3 = gottcha2.parse_taxids('all', mock_taxa_df, "mock_file.tsv")
            
            self.assertEqual(len(qualified_taxids3), 2)
            self.assertIn('12345', qualified_taxids3)
            self.assertIn('67890', qualified_taxids3)

if __name__ == '__main__':
    unittest.main()