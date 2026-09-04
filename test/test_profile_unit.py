import io
import os
import sys
import tempfile
import types
import unittest
from contextlib import redirect_stderr


if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")

from gottcha.utils import profile


class TestProfileUtils(unittest.TestCase):
    def test_parse_args_defaults_for_short_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix])

            self.assertEqual(args.dbLevel, "species")
            self.assertEqual(args.matchIdentity, 0.95)
            self.assertEqual(args.errorRate, 0.005)
            self.assertEqual(args.prefix, "reads")
            self.assertEqual(args.input[0], os.path.abspath(read_path))

    def test_parse_args_nanopore_defaults_to_direct_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "ont.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix, "-np"])

            self.assertTrue(args.nanopore)
            self.assertFalse(args.ont_chunk)
            self.assertEqual(args.matchIdentity, 0.85)
            self.assertEqual(args.matchFraction, 0.05)
            self.assertEqual(args.matchLength, 100)
            self.assertEqual(args.errorRate, 0.03)
            self.assertEqual(args.presetx, "lr:hq")
            self.assertEqual(args.m2options, "-n1 -m25 -s100 --no-long-join")
            self.assertEqual(args.ont_max_secondary, 30)
            self.assertEqual(args.ont_secondary_ratio, 0.5)
            self.assertEqual(args.ont_min_species_support, 0.6)

    def test_parse_args_nanopore_chunk_mode_uses_chunk_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "ont.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args(
                "test",
                ["profile", "-i", read_path, "-d", db_prefix, "-np", "--ont-chunk"],
            )

            self.assertTrue(args.nanopore)
            self.assertTrue(args.ont_chunk)
            self.assertEqual(args.matchIdentity, 0.85)
            self.assertEqual(args.matchFraction, 0.85)
            self.assertEqual(args.matchLength, 100)
            self.assertEqual(args.errorRate, 0.03)
            self.assertEqual(args.presetx, "sr")
            self.assertEqual(args.m2options, "-s120")

    def test_parse_args_rejects_ont_chunk_without_nanopore_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            stderr = io.StringIO()
            with redirect_stderr(stderr), self.assertRaises(SystemExit):
                profile.parse_args(
                    "test",
                    ["profile", "-i", read_path, "-d", db_prefix, "--ont-chunk"],
                )

            self.assertIn("--ont-chunk requires --nanopore", stderr.getvalue())

    def test_parse_args_extractfullref_and_nocutoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args(
                "test",
                ["profile", "-i", read_path, "-d", db_prefix, "-ef", "-nc"],
            )
            self.assertEqual(args.extract, "all:20:fasta")
            self.assertEqual(args.sniScore, "0,0,0")

    def test_load_database_stats_with_header_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats_path = os.path.join(tmp, "db.stats")
            with open(stats_path, "w") as f:
                f.write("Rank\tName\tTaxid\tSK\tNum\tMax\tMin\tTotalLength\tGenomeSize\tNote\n")
                f.write("species\tX\t123\tB\t1\t0\t0\t1000\t1200\tok\n")

            df = profile.load_database_stats(stats_path)
            self.assertIn("123", df.index)
            self.assertEqual(int(df.loc["123", "TotalLength"]), 1000)
            self.assertEqual(int(df.loc["123", "GenomeSize"]), 1200)

    def test_load_acc_list_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "acc.txt")
            open(p, "w").close()
            self.assertEqual(profile.load_acc_list(p), set())


if __name__ == "__main__":
    unittest.main()
