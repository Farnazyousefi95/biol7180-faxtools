"""
Unit tests for faxtools.

Run all tests from the repository root with::

    python3 -m unittest -v

Each test method is named to describe the behaviour it verifies, which makes
the verbose output readable.
"""

import os
import sys
import tempfile
import unittest

# ---------------------------------------------------------------------------
# Make sure the repo root is on sys.path so we can import faxtools directly
# regardless of the working directory.
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import faxtools  # noqa: E402 — import after path setup


# ---------------------------------------------------------------------------
# Helper: create temporary FASTA files for tests
# ---------------------------------------------------------------------------

def _write_temp_fasta(content):
    """Write *content* to a temporary .fasta file and return the path.

    The caller is responsible for deleting the file when done (or using
    ``addCleanup`` in a test case).
    """
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    )
    tmp.write(content)
    tmp.close()
    return tmp.name


# ===================================================================
# Tests for parse_fasta
# ===================================================================

class TestParseFasta(unittest.TestCase):
    """Tests for the ``parse_fasta`` function."""

    def test_single_sequence(self):
        """A file with one header and one sequence line."""
        path = _write_temp_fasta(">seq1 description\nACGTACGT\n")
        self.addCleanup(os.remove, path)

        records = list(faxtools.parse_fasta(path))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0], ("seq1 description", "ACGTACGT"))

    def test_multiple_sequences(self):
        """A file with several records."""
        content = ">s1\nAAAA\n>s2\nCCCC\n>s3\nGGGG\n"
        path = _write_temp_fasta(content)
        self.addCleanup(os.remove, path)

        records = list(faxtools.parse_fasta(path))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[1][1], "CCCC")

    def test_multiline_sequence(self):
        """Sequences split across multiple lines are concatenated."""
        content = ">s1\nACGT\nTGCA\nAAAA\n"
        path = _write_temp_fasta(content)
        self.addCleanup(os.remove, path)

        records = list(faxtools.parse_fasta(path))
        self.assertEqual(records[0][1], "ACGTTGCAAAAA")

    def test_blank_lines_ignored(self):
        """Blank lines between records should not cause errors."""
        content = "\n>s1\nACGT\n\n>s2\nTGCA\n\n"
        path = _write_temp_fasta(content)
        self.addCleanup(os.remove, path)

        records = list(faxtools.parse_fasta(path))
        self.assertEqual(len(records), 2)

    def test_missing_file_raises(self):
        """FileNotFoundError if the path does not exist."""
        with self.assertRaises(FileNotFoundError):
            list(faxtools.parse_fasta("/no/such/file.fasta"))

    def test_no_header_raises(self):
        """ValueError if the first non-blank line is not a header."""
        path = _write_temp_fasta("ACGTACGT\n>late_header\nGGGG\n")
        self.addCleanup(os.remove, path)

        with self.assertRaises(ValueError):
            list(faxtools.parse_fasta(path))


# ===================================================================
# Tests for calculate_gc_percent
# ===================================================================

class TestCalculateGcPercent(unittest.TestCase):
    """Tests for the ``calculate_gc_percent`` function."""

    def test_all_gc(self):
        """A sequence of only G and C should be 100%."""
        self.assertAlmostEqual(faxtools.calculate_gc_percent("GCGCGC"), 100.0)

    def test_all_at(self):
        """A sequence of only A and T should be 0%."""
        self.assertAlmostEqual(faxtools.calculate_gc_percent("ATATAT"), 0.0)

    def test_mixed(self):
        """50/50 mix should be 50%."""
        self.assertAlmostEqual(faxtools.calculate_gc_percent("ACGT"), 50.0)

    def test_case_insensitive(self):
        """Lowercase bases should be counted correctly."""
        self.assertAlmostEqual(faxtools.calculate_gc_percent("acgt"), 50.0)

    def test_gaps_ignored(self):
        """Gap characters ('-') should not affect the percentage."""
        self.assertAlmostEqual(
            faxtools.calculate_gc_percent("A-C-G-T"), 50.0
        )

    def test_empty_sequence(self):
        """An empty sequence should return 0.0 without error."""
        self.assertAlmostEqual(faxtools.calculate_gc_percent(""), 0.0)


# ===================================================================
# Tests for sequence_length
# ===================================================================

class TestSequenceLength(unittest.TestCase):
    """Tests for the ``sequence_length`` function."""

    def test_no_gaps(self):
        """Length of a gap-free sequence."""
        self.assertEqual(faxtools.sequence_length("ACGTACGT"), 8)

    def test_with_gaps(self):
        """Gaps should be excluded from the length."""
        self.assertEqual(faxtools.sequence_length("AC--GT"), 4)

    def test_empty(self):
        """Empty string has length 0."""
        self.assertEqual(faxtools.sequence_length(""), 0)


# ===================================================================
# Tests for summarize_fasta
# ===================================================================

class TestSummarizeFasta(unittest.TestCase):
    """Tests for the ``summarize_fasta`` function."""

    def test_basic_summary(self):
        """Verify counts and statistics for a known input."""
        content = ">s1\nAAAA\n>s2\nCCCCCCCC\n>s3\nGG\n"
        path = _write_temp_fasta(content)
        self.addCleanup(os.remove, path)

        result = summarize = faxtools.summarize_fasta(path)

        self.assertEqual(summarize["n_seqs"], 3)
        self.assertEqual(summarize["total_bp"], 14)   # 4 + 8 + 2
        self.assertEqual(summarize["min_len"], 2)
        self.assertEqual(summarize["max_len"], 8)
        self.assertAlmostEqual(summarize["mean_len"], 4.7, places=1)

    def test_gc_percent_in_summary(self):
        """GC% should reflect overall base composition."""
        # 4 A's + 8 C's + 2 G's = 10 GC out of 14 total
        content = ">s1\nAAAA\n>s2\nCCCCCCCC\n>s3\nGG\n"
        path = _write_temp_fasta(content)
        self.addCleanup(os.remove, path)

        result = faxtools.summarize_fasta(path)
        expected_gc = (10 / 14) * 100.0
        self.assertAlmostEqual(result["gc_percent"], round(expected_gc, 2))

    def test_empty_file_raises(self):
        """A FASTA file with no sequences should raise ValueError."""
        path = _write_temp_fasta("")
        self.addCleanup(os.remove, path)

        with self.assertRaises(ValueError):
            faxtools.summarize_fasta(path)


# ===================================================================
# Tests for filter_sequences
# ===================================================================

class TestFilterSequences(unittest.TestCase):
    """Tests for the ``filter_sequences`` function."""

    def setUp(self):
        """Create a temp FASTA with sequences of known lengths."""
        self.content = ">short\nAA\n>medium\nACGTACGT\n>long\nACGTACGTACGTACGT\n"
        self.path = _write_temp_fasta(self.content)
        self.addCleanup(os.remove, self.path)

    def test_filter_keeps_long_sequences(self):
        """Only sequences >= min_len should pass."""
        kept = list(faxtools.filter_sequences([self.path], min_len=8))
        lengths = [faxtools.sequence_length(seq) for _, seq in kept]
        self.assertTrue(all(l >= 8 for l in lengths))
        self.assertEqual(len(kept), 2)  # medium (8) and long (16)

    def test_filter_keeps_all_when_threshold_is_low(self):
        """All sequences should pass when min_len is 1."""
        kept = list(faxtools.filter_sequences([self.path], min_len=1))
        self.assertEqual(len(kept), 3)

    def test_filter_excludes_all_when_threshold_is_high(self):
        """No sequences should pass if min_len exceeds all lengths."""
        kept = list(faxtools.filter_sequences([self.path], min_len=100))
        self.assertEqual(len(kept), 0)


# ===================================================================
# Tests for write_fasta
# ===================================================================

class TestWriteFasta(unittest.TestCase):
    """Tests for the ``write_fasta`` function."""

    def test_roundtrip(self):
        """Writing and re-reading should preserve headers and sequences."""
        records = [("s1", "AAAA"), ("s2", "CCCCGGGG")]
        tmp = tempfile.NamedTemporaryFile(
            suffix=".fasta", delete=False, mode="w"
        )
        tmp.close()
        self.addCleanup(os.remove, tmp.name)

        faxtools.write_fasta(records, tmp.name)
        recovered = list(faxtools.parse_fasta(tmp.name))

        self.assertEqual(len(recovered), 2)
        self.assertEqual(recovered[0], ("s1", "AAAA"))
        self.assertEqual(recovered[1], ("s2", "CCCCGGGG"))


# ===================================================================
# Tests for format_summary_tsv
# ===================================================================

class TestFormatSummaryTsv(unittest.TestCase):
    """Tests for the ``format_summary_tsv`` function."""

    def test_header_and_rows(self):
        """Output should have a header line followed by data rows."""
        summaries = [
            {
                "file": "test.fasta",
                "n_seqs": 3,
                "total_bp": 100,
                "min_len": 10,
                "mean_len": 33.3,
                "max_len": 50,
                "gc_percent": 45.12,
            }
        ]
        tsv = faxtools.format_summary_tsv(summaries)
        lines = tsv.strip().split("\n")

        self.assertEqual(len(lines), 2)  # header + 1 data row
        self.assertTrue(lines[0].startswith("file\t"))
        self.assertIn("test.fasta", lines[1])


# ===================================================================
# Tests for the CLI (main function)
# ===================================================================

class TestCLI(unittest.TestCase):
    """Integration tests that exercise the ``main`` entry point."""

    def setUp(self):
        """Create a small temp FASTA for CLI tests."""
        self.content = ">s1\nACGTACGT\n>s2\nGGGGGGGGGGGGGGGG\n"
        self.path = _write_temp_fasta(self.content)
        self.addCleanup(os.remove, self.path)

    def test_main_stdout(self):
        """main() should print a TSV summary to stdout."""
        # Capture stdout
        from io import StringIO
        captured = StringIO()
        old_stdout = sys.stdout
        sys.stdout = captured
        try:
            faxtools.main([self.path])
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        self.assertIn("file\t", output)
        self.assertIn("2", output)  # n_seqs

    def test_main_out_file(self):
        """--out should write the TSV to a file."""
        out_tmp = tempfile.NamedTemporaryFile(
            suffix=".tsv", delete=False, mode="w"
        )
        out_tmp.close()
        self.addCleanup(os.remove, out_tmp.name)

        faxtools.main([self.path, "--out", out_tmp.name])

        with open(out_tmp.name) as fh:
            content = fh.read()
        self.assertIn("file\t", content)

    def test_main_filter(self):
        """--min-len + --filtered-out should produce a filtered FASTA."""
        filt_tmp = tempfile.NamedTemporaryFile(
            suffix=".fasta", delete=False, mode="w"
        )
        filt_tmp.close()
        self.addCleanup(os.remove, filt_tmp.name)

        faxtools.main([
            self.path, "--min-len", "10", "--filtered-out", filt_tmp.name
        ])

        records = list(faxtools.parse_fasta(filt_tmp.name))
        # Only s2 (length 16) should pass the min-len=10 filter
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][0], "s2")


# ===================================================================
# Tests using the bundled example data files
# ===================================================================

class TestExampleData(unittest.TestCase):
    """Smoke tests that run against the files in data/."""

    DATA_DIR = os.path.join(REPO_ROOT, "data")

    def test_example_seqs1_parses(self):
        """data/example-seqs1.fasta should parse without error."""
        path = os.path.join(self.DATA_DIR, "example-seqs1.fasta")
        if not os.path.isfile(path):
            self.skipTest("example-seqs1.fasta not found")
        records = list(faxtools.parse_fasta(path))
        self.assertGreater(len(records), 0)

    def test_example_seqs2_parses(self):
        """data/example-seqs2.fasta should parse without error."""
        path = os.path.join(self.DATA_DIR, "example-seqs2.fasta")
        if not os.path.isfile(path):
            self.skipTest("example-seqs2.fasta not found")
        records = list(faxtools.parse_fasta(path))
        self.assertGreater(len(records), 0)

    def test_example_seqs1_summary(self):
        """Summarize data/example-seqs1.fasta and check basic fields."""
        path = os.path.join(self.DATA_DIR, "example-seqs1.fasta")
        if not os.path.isfile(path):
            self.skipTest("example-seqs1.fasta not found")
        summary = faxtools.summarize_fasta(path)
        self.assertEqual(summary["n_seqs"], 5)
        self.assertGreater(summary["total_bp"], 0)
        self.assertGreater(summary["gc_percent"], 0)


if __name__ == "__main__":
    unittest.main()
