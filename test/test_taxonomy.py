import unittest
import os
import gottcha.scripts.taxonomy as taxonomy

class TestTaxonomy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_data_dir = os.path.dirname(os.path.abspath(__file__))
        taxonomy_dir = os.path.join(test_data_dir, "gottcha2_database_test/taxonomy_db")
        taxonomy.loadTaxonomy(taxonomy_dir)

    def test_taxid2rank(self):
        rank = taxonomy.taxid2rank("1")
        self.assertEqual(rank, "root")

    def test_taxid2name(self):
        name = taxonomy.taxid2name("1")
        self.assertEqual(name, "root")

    def test_taxid2parent(self):
        parent = taxonomy.taxid2parent("1")
        self.assertEqual(parent, "1")

    def test_is_descendant(self):
        result = taxonomy._is_descendant("1", "1")
        self.assertFalse(result)

    def test_taxid2lineage(self):
        lineage = taxonomy.taxid2lineage("1")
        self.assertIsInstance(lineage, str)

if __name__ == '__main__':
    unittest.main()
