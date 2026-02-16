# biol7180-faxtools

**BIOL 7180 (Scripting for Biologists) — Solo Project**
**Student:** Farnaz Yousefi

## Project overview

`faxtools` is a lightweight command-line tool to **summarize FASTA files** and (optionally) **filter sequences by minimum length**. It is designed for quick QC of biological sequence datasets (e.g., amplicons, reference genes, contigs).

## Features (MVP)

* Robust FASTA parsing (supports multi-line sequences)
* Per-file summary statistics:

  * number of sequences
  * total bases
  * min / mean / max sequence length
  * GC% (case-insensitive; ignores gaps `-`)
* Optional filtering:

  * keep sequences with length ≥ `--min-len`
  * write filtered sequences to a new FASTA file

## Repository contents

* `faxtools.py` — main script
* `data/` — example FASTA files for testing/demo
* `tests/` — automated tests
* `PROPOSAL.md` — project proposal (due Feb 18)

## Requirements

* Python 3.8+ (no external packages required)

## Usage (examples)

Summarize one FASTA:

```bash
python3 faxtools.py data/example-seqs1.fasta
```

Summarize multiple FASTA files and save a TSV:

```bash
python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta --out summary.tsv
```

Filter sequences by minimum length and write filtered FASTA:

```bash
python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta --min-len 90 --filtered-out filtered.fasta
```

## Run tests

```bash
python3 -m unittest -v
```

## Expected output format

The summary output is tab-separated with a header:

```
file	n_seqs	total_bp	min_len	mean_len	max_len	gc_percent
```

## GitHub link

https://github.com/Farnazyousefi95/biol7180-faxtools

