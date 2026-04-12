# faxtools — FASTA Summary & Filtering Tool

**BIOL 7180 (Scripting for Biologists) — Solo Project**  
**Student:** Farnaz Yousefi  
**Semester:** Spring 2026  
**GitHub repo:** <https://github.com/Farnazyousefi95/biol7180-faxtools>

---

## Table of Contents

1. [Project Overview](#project-overview)  
2. [Biological Motivation](#biological-motivation)  
3. [Features](#features)  
4. [Repository Contents](#repository-contents)  
5. [Requirements](#requirements)  
6. [Installation](#installation)  
7. [Tutorial](#tutorial)  
   - [1 — Summarize a single FASTA file](#1--summarize-a-single-fasta-file)  
   - [2 — Summarize multiple files](#2--summarize-multiple-files)  
   - [3 — Save the summary as a TSV](#3--save-the-summary-as-a-tsv)  
   - [4 — Filter sequences by minimum length](#4--filter-sequences-by-minimum-length)  
   - [5 — Combine summary + filtering in one command](#5--combine-summary--filtering-in-one-command)  
8. [Interpreting the Output](#interpreting-the-output)  
9. [Running the Tests](#running-the-tests)  
10. [Code Design and Best Practices](#code-design-and-best-practices)  
11. [Limitations and Future Work](#limitations-and-future-work)  
12. [License](#license)

---

## Project Overview

`faxtools` is a lightweight, zero-dependency, command-line Python tool that
summarizes FASTA files and optionally filters sequences by minimum length. It
is designed for **quick quality control (QC)** of biological sequence datasets
such as amplicons, reference genes, or assembled contigs.

Given one or more FASTA files, `faxtools` reports:

| Metric | Description |
|---|---|
| `n_seqs` | Number of sequences in the file |
| `total_bp` | Total bases (gaps excluded) |
| `min_len` | Shortest sequence length |
| `mean_len` | Mean sequence length |
| `max_len` | Longest sequence length |
| `gc_percent` | Overall GC content (%) |

## Biological Motivation

FASTA is the standard format for nucleotide and protein sequences. Before
running downstream analyses (e.g., alignment, phylogenetics, assembly
evaluation), researchers commonly need to answer quick QC questions:

- How many sequences are in this file?
- What is the length distribution?
- What is the GC content? (useful for detecting contamination or confirming
  species identity)
- Are there very short sequences that should be filtered out?

`faxtools` answers all of these questions with a single command and produces
tab-separated output that is easy to paste into a spreadsheet or pipe into
other tools.

## Features

- **Robust FASTA parsing** — handles multi-line sequences, blank lines, and
  validates that each file begins with a `>` header.
- **Per-file summary statistics** — sequence count, total bases, min/mean/max
  length, and GC%.
- **Minimum-length filtering** — keep only sequences ≥ a threshold and write
  them to a new FASTA file.
- **TSV output** — summary table to stdout or to a file (`--out`).
- **No external dependencies** — pure Python 3.8+ standard library.

## Repository Contents

```
biol7180-faxtools/
├── faxtools.py                 # Main script
├── data/
│   ├── example-seqs1.fasta     # 5 sequences (diverse organisms)
│   └── example-seqs2.fasta     # 5 contigs (including gaps & GC-rich)
├── tests/
│   ├── __init__.py
│   └── test_faxtools.py        # Unit & integration tests
├── PROPOSAL.md                 # Original project proposal
└── README.md                   # This file (tutorial + documentation)
```

## Requirements

- **Python 3.8** or newer (tested on 3.10 and 3.12)
- No external packages — only the Python standard library is used.

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/Farnazyousefi95/biol7180-faxtools.git
   cd biol7180-faxtools
   ```

2. **Verify Python version:**

   ```bash
   python3 --version   # should be 3.8+
   ```

3. **That's it!** There is nothing to install. The script runs directly.

---

## Tutorial

The examples below use the two FASTA files bundled in the `data/` directory.
You can follow along by running each command from the repository root.

### 1 — Summarize a single FASTA file

```bash
python3 faxtools.py data/example-seqs1.fasta
```

**Expected output** (tab-separated; columns are aligned here for readability):

```
file               n_seqs  total_bp  min_len  mean_len  max_len  gc_percent
example-seqs1.fasta  5       1019      30       203.8     430      51.52
```

> **What this tells you:** The file contains 5 sequences ranging from 30 bp to
> 430 bp, with an overall GC content of about 52%.

### 2 — Summarize multiple files

```bash
python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta
```

This prints one row per file, making it easy to compare datasets at a glance.

### 3 — Save the summary as a TSV

```bash
python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta --out summary.tsv
```

The file `summary.tsv` is now ready to open in Excel, R, or any
spreadsheet program. A confirmation message is printed to stderr:

```
Summary written to summary.tsv
```

### 4 — Filter sequences by minimum length

Suppose you want to discard sequences shorter than 90 bp:

```bash
python3 faxtools.py data/example-seqs1.fasta --min-len 90 --filtered-out filtered.fasta
```

This will:
1. Print the summary table to stdout (as usual).
2. Write only the sequences with length ≥ 90 to `filtered.fasta`.
3. Report how many sequences passed the filter (to stderr).

You can verify the filtered file:

```bash
grep -c "^>" filtered.fasta
```

### 5 — Combine summary + filtering in one command

```bash
python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta \
    --out summary.tsv \
    --min-len 90 \
    --filtered-out filtered.fasta
```

This saves the summary **and** writes the filtered FASTA in a single run.

---

## Interpreting the Output

| Column | Meaning | Why it matters |
|---|---|---|
| `file` | Input file name | Identifies the source dataset |
| `n_seqs` | Number of sequences | Confirms expected count after demultiplexing or assembly |
| `total_bp` | Total non-gap bases | Useful for estimating genome coverage |
| `min_len` | Shortest sequence | Flags unexpectedly short reads or contigs |
| `mean_len` | Average length | Quick summary of length distribution |
| `max_len` | Longest sequence | Sanity check on assembly or amplicon size |
| `gc_percent` | Overall GC% | Can flag contamination if GC deviates from expected range |

---

## Running the Tests

From the repository root, run:

```bash
python3 -m unittest -v
```

You should see output similar to:

```
test_single_sequence (tests.test_faxtools.TestParseFasta) ... ok
test_multiple_sequences (tests.test_faxtools.TestParseFasta) ... ok
test_multiline_sequence (tests.test_faxtools.TestParseFasta) ... ok
...
----------------------------------------------------------------------
Ran 24 tests in 0.02s

OK
```

The test suite covers:

- **FASTA parsing:** single/multi-sequence files, multi-line sequences, blank
  lines, missing files, and invalid formats.
- **GC calculation:** pure GC, pure AT, mixed, case-insensitive, gaps, and
  empty sequences.
- **Sequence length:** with and without gaps.
- **Summary statistics:** counts, lengths, GC% for known inputs.
- **Filtering:** threshold boundary cases (all pass, some pass, none pass).
- **FASTA writing:** round-trip (write then re-read) integrity.
- **TSV formatting:** header presence and data content.
- **CLI integration:** stdout output, `--out` file writing, and
  `--min-len`/`--filtered-out` filtering.
- **Example data smoke tests:** ensures the bundled data files parse and
  summarize without errors.

---

## Code Design and Best Practices

The code was written to follow the best practices taught in BIOL 7180:

- **Modularity:** Each task (parsing, GC calculation, summarizing, filtering,
  writing, TSV formatting) is in its own function. Functions are small and
  do one thing.
- **Expressive naming:** Variable and function names describe their purpose
  (e.g., `calculate_gc_percent`, `sequence_length`, `filter_sequences`).
- **Docstrings:** Every function has a NumPy-style docstring documenting its
  parameters, return values, and raised exceptions.
- **Error handling:** The parser validates that files exist and that headers
  appear before sequence data. The summarizer raises a clear error on empty
  files.
- **Testability:** `main()` accepts an `argv` parameter so the CLI can be
  tested programmatically without subprocess calls.
- **No external dependencies:** The script uses only the Python standard
  library, making it easy to run anywhere Python 3.8+ is available.

---

## Limitations and Future Work

- Currently only handles DNA sequences (does not distinguish RNA or protein).
- N50 / N90 statistics could be added for assembly QC.
- A `--format` option (e.g., CSV, JSON) could provide more flexible output.
- Integration with BioPython could enable additional format support (GenBank,
  FASTQ).

---

## License

This project was developed for educational purposes as part of BIOL 7180 at
Auburn University. Feel free to use and modify it for your own work.

---

## Getting Help

```bash
python3 faxtools.py --help
```

```
usage: faxtools.py [-h] [--out TSV] [--min-len N] [--filtered-out FASTA]
                   FASTA [FASTA ...]

Summarize one or more FASTA files (sequence counts, length statistics, GC%)
and optionally filter sequences by minimum length.

positional arguments:
  FASTA                 One or more FASTA files to process.

options:
  -h, --help            show this help message and exit
  --out TSV             Write the summary table to this TSV file (default: stdout).
  --min-len N           Minimum sequence length for filtering (requires --filtered-out).
  --filtered-out FASTA  Write sequences passing --min-len to this FASTA file.

examples:
  python3 faxtools.py data/example-seqs1.fasta
  python3 faxtools.py *.fasta --out summary.tsv
  python3 faxtools.py *.fasta --min-len 100 --filtered-out long_seqs.fasta
```
