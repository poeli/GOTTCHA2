import unittest
import os
import sys
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from pathlib import Path
from io import StringIO

# Get the path to the project root directory
current_dir = Path(__file__).parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import the module
import src.gottcha2.gottcha2 as gottcha2

class TestGottcha2CLI(unittest.TestCase):
    """Test command-line interface functionality of GOTTCHA2."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        # Mock database file
        self.db_path = os.path.join(self.test_dir, "test_db")
        # Create a dummy .mmi file to pass database existence check
        with open(f"{self.db_path}.mmi", 'w') as f:
            f.write("dummy minimap2 index")
        # Create a dummy stats file
        with open(f"{self.db_path}.stats", 'w') as f:
            f.write("species\tE. coli\t12345\tBacteria\t1\t100\t50\t1000\n")
    
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
        test_fastq = os.path.join(self.test_dir, "test.fastq")
        with open(test_fastq, 'w') as f:
            f.write("@seq1\nACGT\n+\nHHHH\n")
        
        args = gottcha2.parse_params("1.0.0", ['-i', test_fastq, '-d', self.db_path, '-l', 'species'])
        self.assertEqual(args.database, self.db_path)
        self.assertEqual(args.dbLevel, 'species')
        
        # Test noCutoff option
        args = gottcha2.parse_params("1.0.0", ['-i', test_fastq, '-d', self.db_path, '-l', 'species', '-nc'])
        self.assertEqual(args.minCov, 0)
        self.assertEqual(args.minReads, 0)
        self.assertEqual(args.minLen, 0)
        self.assertEqual(args.matchFactor, 0)
        self.assertEqual(args.maxZscore, 0)
        
        # Test nanopore option
        args = gottcha2.parse_params("1.0.0", ['-i', test_fastq, '-d', self.db_path, '-l', 'species', '-np'])
        self.assertEqual(args.presetx, 'map-ont')
        self.assertEqual(args.minReads, 1)
        self.assertEqual(args.matchFactor, 0)
    
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
    
    @patch('src.gottcha2.gottcha2.subprocess.Popen')
    def test_readMapping(self, mock_popen):
        """Test read mapping function."""
        # Mock the subprocess.Popen
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("", "")
        mock_process.poll.return_value = 0  # Success exit code
        mock_popen.return_value = mock_process
        
        # Create test files
        test_fastq = os.path.join(self.test_dir, "test.fastq")
        with open(test_fastq, 'w') as f:
            f.write("@seq1\nACGT\n+\nHHHH\n")
        
        # Create a file-like object for input reads
        class MockFile:
            def __init__(self, name):
                self.name = name
        
        reads = [MockFile(test_fastq)]
        logfile = os.path.join(self.test_dir, "test.log")
        samfile = os.path.join(self.test_dir, "test.sam")
        
        # Test the function
        exit_code, cmd, errs = gottcha2.readMapping(
            reads, self.db_path, 1, 10, 'sr', samfile, logfile, False
        )
        
        # Check the results
        self.assertEqual(exit_code, 0)
        self.assertIn("minimap2", cmd)
        self.assertIn(self.db_path, cmd)
        mock_popen.assert_called_once()

if __name__ == '__main__':
    unittest.main()
