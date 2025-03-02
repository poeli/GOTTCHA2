import unittest
import os
import sys
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from pathlib import Path
import json

# Get the path to the project root directory
current_dir = Path(__file__).parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import the module - adjust if needed based on project structure
try:
    from gottcha.scripts import taxonomy as gt
except ImportError:
    import gottcha.scripts.taxonomy as gt

class TestTaxonomy(unittest.TestCase):
    """Tests for taxonomy.py functions."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment once for all tests."""
        # Create a temporary directory for test data
        cls.test_dir = tempfile.mkdtemp()
        
        # Create minimal taxonomy files for testing
        cls.create_test_taxonomy_files()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests."""
        # Remove temporary directory and files
        shutil.rmtree(cls.test_dir)
    
    @classmethod
    def create_test_taxonomy_files(cls):
        """Create minimalist taxonomy files for testing."""
        # Create the taxonomy directory
        taxonomy_dir = os.path.join(cls.test_dir, "taxonomy_db")
        os.makedirs(taxonomy_dir, exist_ok=True)
        
        # Create names.dmp
        with open(os.path.join(taxonomy_dir, "names.dmp"), 'w') as f:
            f.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
            f.write("2\t|\tBacteria\t|\t\t|\tscientific name\t|\n")
            f.write("6\t|\tProteobacteria\t|\t\t|\tscientific name\t|\n")
            f.write("1224\t|\tProteobacteria\t|\t\t|\tscientific name\t|\n")  # Duplicate name
            f.write("1236\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|\n")
            f.write("91347\t|\tEnterobacterales\t|\t\t|\tscientific name\t|\n")
            f.write("543\t|\tEnterobacteriaceae\t|\t\t|\tscientific name\t|\n")
            f.write("561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n")
            f.write("562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n")
            f.write("11111\t|\tE. coli strain K12\t|\t\t|\tscientific name\t|\n")
        
        # Create nodes.dmp
        with open(os.path.join(taxonomy_dir, "nodes.dmp"), 'w') as f:
            f.write("1\t|\t1\t|\troot\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # root
            f.write("2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Bacteria
            f.write("6\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Proteobacteria old
            f.write("1224\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Proteobacteria
            f.write("1236\t|\t1224\t|\tclass\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Gammaproteobacteria
            f.write("91347\t|\t1236\t|\torder\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Enterobacterales
            f.write("543\t|\t91347\t|\tfamily\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Enterobacteriaceae
            f.write("561\t|\t543\t|\tgenus\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Escherichia
            f.write("562\t|\t561\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # Escherichia coli
            f.write("11111\t|\t562\t|\tno rank\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")  # E. coli strain
        
        # Create merged.dmp
        with open(os.path.join(taxonomy_dir, "merged.dmp"), 'w') as f:
            f.write("6\t|\t1224\t|\n")  # Old Proteobacteria merged into new Proteobacteria
        
        # Create major_level_to_abbr.json
        with open(os.path.join(taxonomy_dir, "major_level_to_abbr.json"), 'w') as f:
            json.dump({
                "superkingdom": "sk",
                "phylum": "p",
                "class": "c",
                "order": "o",
                "family": "f",
                "genus": "g",
                "species": "s",
                "strain": "n"
            }, f)
    
    def setUp(self):
        """Set up test environment for each test."""
        # Reset taxonomy globals before each test to ensure clean state
        gt.taxDepths = {"1": 0}
        gt.taxParents = {"1": "1"}
        gt.taxRanks = {"1": "root"}
        gt.taxNames = {"1": "root"}
        gt.taxMerged = {}
        gt.taxNumChilds = {}
        gt.accTid = {}
        gt.tidLineage = {}
        gt.tidLineageDict = {}
        gt.nameTid = {}
        gt.major_level_to_abbr = {}
        gt.abbr_to_major_level = {}
        
        # Load the test taxonomy
        taxonomy_dir = os.path.join(self.test_dir, "taxonomy_db")
        gt.loadTaxonomy(dbpath=taxonomy_dir, auto_download=False)
    
    def test_taxonomy_loading(self):
        """Test that taxonomy data is correctly loaded."""
        # Verify some key items were loaded
        self.assertEqual(gt.taxNames["1"], "root")
        self.assertEqual(gt.taxNames["562"], "Escherichia coli")
        self.assertEqual(gt.taxRanks["562"], "species")
        self.assertEqual(gt.taxRanks["1"], "root")
        self.assertEqual(gt.taxParents["562"], "561")  # E. coli's parent should be Escherichia
        
        # Check if merged taxids were loaded
        self.assertEqual(gt.taxMerged["6"], "1224")  # Old Proteobacteria merged into new Proteobacteria
    
    def test_taxid2rank(self):
        """Test getting rank from taxid."""
        self.assertEqual(gt.taxid2rank("1"), "root")
        self.assertEqual(gt.taxid2rank("2"), "superkingdom")
        self.assertEqual(gt.taxid2rank("1224"), "phylum")
        self.assertEqual(gt.taxid2rank("562"), "species")
        self.assertEqual(gt.taxid2rank("11111"), "strain")  # Test a no rank that should be guessed as strain
        
        # Test merged taxid
        self.assertEqual(gt.taxid2rank("6"), "phylum")  # Should get rank from the merged target (1224)
        
        # Test invalid taxid
        self.assertEqual(gt.taxid2rank("999999"), "unknown")
    
    def test_taxid2name(self):
        """Test getting name from taxid."""
        self.assertEqual(gt.taxid2name("1"), "root")
        self.assertEqual(gt.taxid2name("562"), "Escherichia coli")
        
        # Test merged taxid
        self.assertEqual(gt.taxid2name("6"), "Proteobacteria")  # Should get name from merged target (1224)
        
        # Test invalid taxid
        self.assertEqual(gt.taxid2name("999999"), "unknown")
    
    def test_taxid2parent(self):
        """Test getting parent taxid."""
        self.assertEqual(gt.taxid2parent("562"), "561")  # E. coli -> Escherichia
        self.assertEqual(gt.taxid2parent("561"), "543")  # Escherichia -> Enterobacteriaceae
        
        # Test with norank=True to get immediate parent
        self.assertEqual(gt.taxid2parent("11111", norank=True), "562")  # E. coli strain -> E. coli
        
        # Test with norank=False (default) to skip no-rank parents
        self.assertEqual(gt.taxid2parent("11111"), "562")  # E. coli strain -> Escherichia (skipping E. coli which is no-rank)
        
        # Test invalid taxid
        self.assertEqual(gt.taxid2parent("999999"), "unknown")
    
    def test_taxid2lineage(self):
        """Test generating lineage string."""
        # Test with default parameters
        lineage = gt.taxid2lineage("562")  # E. coli
        self.assertIn("superkingdom|2|Bacteria", lineage)
        self.assertIn("phylum|1224|Proteobacteria", lineage)
        self.assertIn("species|562|Escherichia coli", lineage)
        
        # Test with space2underscore=True
        lineage = gt.taxid2lineage("562", space2underscore=True)
        self.assertIn("species|562|Escherichia_coli", lineage)
        
        # Test with custom separator
        lineage = gt.taxid2lineage("562", sep=";")
        self.assertIn("sk__Bacteria", lineage)
        self.assertIn("s__Escherichia coli", lineage)
        
        # Test invalid taxid
        lineage = gt.taxid2lineage("999999")
        self.assertEqual(lineage, "")
    
    def test_taxidIsLeaf(self):
        """Test checking if a taxid is a leaf node."""
        self.assertTrue(gt.taxidIsLeaf("11111"))  # E. coli strain is a leaf
        self.assertFalse(gt.taxidIsLeaf("562"))   # E. coli is not a leaf (has strain child)
        
        # Test invalid taxid
        self.assertFalse(gt.taxidIsLeaf("999999"))
    
    def test_taxid2nameOnRank(self):
        """Test getting name at specific rank."""
        # Get name at different ranks
        self.assertEqual(gt.taxid2nameOnRank("11111", "species"), "Escherichia coli")
        self.assertEqual(gt.taxid2nameOnRank("11111", "genus"), "Escherichia")
        self.assertEqual(gt.taxid2nameOnRank("11111", "family"), "Enterobacteriaceae")
        
        # Test target rank not in lineage
        self.assertEqual(gt.taxid2nameOnRank("11111", "subphylum"), None)
        
        # Test invalid taxid
        self.assertEqual(gt.taxid2nameOnRank("999999", "species"), "unknown")
    
    def test_taxid2taxidOnRank(self):
        """Test getting taxid at specific rank."""
        # Get taxids at different ranks
        self.assertEqual(gt.taxid2taxidOnRank("11111", "species"), "562")
        self.assertEqual(gt.taxid2taxidOnRank("11111", "genus"), "561")
        
        # Test target rank not in lineage
        self.assertEqual(gt.taxid2taxidOnRank("11111", "subphylum"), None)
        
        # Test invalid taxid
        self.assertEqual(gt.taxid2taxidOnRank("999999", "species"), "unknown")
    
    @patch('builtins.print')  # To suppress print statements
    def test_acc2taxid(self, mock_print):
        """Test converting accession to taxid."""
        # Set up a custom mock file for acc2taxid testing
        acc2taxid_dir = os.path.join(self.test_dir, "taxonomy_db", "accession2taxid")
        os.makedirs(acc2taxid_dir, exist_ok=True)
        
        with open(os.path.join(acc2taxid_dir, "nucl_gb.accession2taxid"), 'w') as f:
            f.write("accession\taccession.version\ttaxid\tgi\n")
            f.write("NC_000913\tNC_000913.3\t562\t0\n")  # E. coli
        
        with patch('gottcha.scripts.taxonomy.acc2taxid_raw', return_value="562"):
            result = gt.acc2taxid("NC_000913")
            self.assertEqual(result, "562")
        
        # Test non-existent accession
        with patch('gottcha.scripts.taxonomy.acc2taxid_raw', return_value=""):
            result = gt.acc2taxid("NONEXISTENT")
            self.assertEqual(result, "")
    
    def test_lca_taxid(self):
        """Test finding lowest common ancestor."""
        # Test LCA of E. coli strains (should be species)
        lca = gt.lca_taxid(["11111", "562"])
        self.assertEqual(lca, "562")  # Should be E. coli
        
        # Test LCA of different species in same genus
        with patch('gottcha.scripts.taxonomy.taxid2lineageDICT', side_effect=[
            # Mocked lineage for 562 (E. coli)
            {"species": {"taxid": "562"}, "genus": {"taxid": "561"}},
            # Mocked lineage for different species in same genus
            {"species": {"taxid": "620"}, "genus": {"taxid": "561"}}
        ]):
            lca = gt.lca_taxid(["562", "620"])
            self.assertEqual(lca, "561")  # Should be genus Escherichia
    
    def test_taxid2fullLineage(self):
        """Test generating full lineage string."""
        lineage = gt.taxid2fullLineage("562")
        
        # Check that the lineage contains all expected levels
        self.assertIn("superkingdom|2|Bacteria", lineage)
        self.assertIn("phylum|1224|Proteobacteria", lineage)
        self.assertIn("species|562|Escherichia_coli", lineage)
        
        # Test with different separator
        lineage = gt.taxid2fullLineage("562", sep=";")
        self.assertIn("phylum__Proteobacteria", lineage)
        
        # Test with abbreviated ranks
        lineage = gt.taxid2fullLineage("562", use_rank_abbr=True)
        self.assertIn("p|1224|Proteobacteria", lineage)

if __name__ == '__main__':
    unittest.main()