# Everything in My Code, Explained

---

## PART 1: MODULES (the import statements at the top)

Every Python script starts by importing modules. These are like toolboxes
that gives us pre-built functionality.

### import argparse

**What:** A module for building command-line interfaces.

**Why you use it:** So users can type things like
`python3 faxtools.py --min-len 90 --out results.tsv` and your script
understands what they want. Without argparse, you'd have to manually
parse the words the user typed.

**Where in your code:** `build_argument_parser()` and `main()`.

### import os

**What:** A module for interacting with the operating system (files,
paths, directories).

**Why you use it:** Two specific functions:
- `os.path.isfile(filepath)` — checks if a file exists before trying
  to open it. Returns True or False.
- `os.path.basename(filepath)` — extracts just the filename from a
  full path. Example: `os.path.basename("/home/user/data/seqs.fasta")`
  returns `"seqs.fasta"`.

**Where in your code:** `parse_fasta()` and `summarize_fasta()`.

### import sys

**What:** A module for interacting with the Python system itself.

**Why you use it:** Two specific things:
- `sys.stderr` — the error output channel (explained earlier).
  Used in `print(f"Error: {err}", file=sys.stderr)`.
- `sys.exit(1)` — stops the script immediately with an error code.
  Code 0 means "everything was fine." Code 1 means "something went
  wrong." Other programs can check this code.

**Where in your code:** `main()` when an error occurs.

### from statistics import mean

**What:** Imports just the `mean` function from Python's statistics
module.

**Why you use it:** To calculate average sequence length. `mean([4, 8, 2])`
returns `4.666...`. You could do `sum(lengths) / len(lengths)` yourself,
but `mean()` is cleaner and handles edge cases.

**Why `from ... import` instead of `import`?** Writing
`from statistics import mean` lets you call `mean()` directly. If you
wrote `import statistics`, you'd need `statistics.mean()` every time.
Both work; this is just shorter.

**Where in your code:** `summarize_fasta()`.

---

## PART 2: EVERY FUNCTION, LINE BY LINE

### Function 1: parse_fasta(filepath)

```python
def parse_fasta(filepath):
```
**`def`** — defines a new function.
**`filepath`** — the parameter (input). When you call
`parse_fasta("data/seqs.fasta")`, filepath becomes `"data/seqs.fasta"`.

```python
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
```
**What:** Checks if the file exists BEFORE trying to open it.
**`raise`** — creates an error on purpose with a helpful message.
**`f"...{filepath}"`** — an f-string. Python replaces `{filepath}` with
the actual value. So if filepath is `"bad.fasta"`, the message becomes
`"FASTA file not found: bad.fasta"`.

```python
    header = None
    sequence_parts = []
```
**What:** Initializes two variables.
- `header = None` — no sequence has been found yet. `None` is Python's
  way of saying "nothing."
- `sequence_parts = []` — an empty list that will collect sequence lines.

```python
    with open(filepath, "r") as fasta_file:
```
**`with`** — a context manager. It opens the file AND automatically
closes it when done (even if an error happens). Without `with`, you'd
need to manually call `fasta_file.close()`.
**`open(filepath, "r")`** — opens the file for reading. `"r"` means
read mode.
**`as fasta_file`** — gives the opened file a name you can use.

```python
        for line_number, raw_line in enumerate(fasta_file, start=1):
```
**`for ... in fasta_file`** — reads the file one line at a time.
**`enumerate(..., start=1)`** — adds a counter. Instead of just getting
each line, you get (1, first_line), (2, second_line), etc. The `start=1`
means counting begins at 1 (not 0). You use this for error messages so
you can say "line 5 has a problem."

```python
            line = raw_line.strip()
```
**`.strip()`** — removes whitespace from both ends of the string.
When Python reads a file, each line ends with `\n` (newline character).
`.strip()` removes that. Example: `"  ACGT\n  ".strip()` becomes `"ACGT"`.

```python
            if not line:
                continue
```
**What:** If the line is empty (blank line), skip it and go to the next
line. `continue` means "jump back to the top of the for loop."

```python
            if line.startswith(">"):
```
**`.startswith(">")`** — returns True if the string begins with `>`.
This is how you detect a FASTA header line.

```python
                if header is not None:
                    yield (header, "".join(sequence_parts))
```
**What:** If we already have a previous sequence stored, hand it back
before starting the new one.
**`yield`** — like `return`, but the function pauses instead of stopping.
Next time you ask for a value, it continues from where it left off.
This makes parse_fasta a **generator**.
**`"".join(sequence_parts)`** — takes a list like
`["ACGT", "TGCA", "AAAA"]` and joins them into one string:
`"ACGTTGCAAAAA"`. The `""` means "put nothing between them."

```python
                header = line[1:].strip()
                sequence_parts = []
```
**`line[1:]`** — slicing. Takes everything AFTER the first character.
So `">seq1 human"` becomes `"seq1 human"` (removes the `>`).
**`sequence_parts = []`** — resets the list for the new sequence.

```python
            else:
                if header is None:
                    raise ValueError(
                        f"{filepath}, line {line_number}: expected a '>' "
                        f"header before sequence data"
                    )
                sequence_parts.append(line)
```
**What:** If the line is NOT a header (no `>`), it must be sequence data.
But if we haven't seen a header yet (`header is None`), the file is
invalid — sequence data appeared before any `>` line. We raise an error.
Otherwise, add this line to our collection.
**`.append(line)`** — adds an item to the end of a list.

```python
        if header is not None:
            yield (header, "".join(sequence_parts))
```
**What:** After the loop ends (end of file), yield the last sequence.
Without this, the final sequence in the file would be lost.

---

### Function 2: calculate_gc_percent(sequence)

```python
    bases = sequence.upper().replace("-", "")
```
**`.upper()`** — converts to uppercase. `"acgt"` becomes `"ACGT"`.
This makes the function case-insensitive.
**`.replace("-", "")`** — removes all gap characters. `"A-C-G-T"`
becomes `"ACGT"`. Method chaining: you can call multiple methods
in a row, left to right.

```python
    if not bases:
        return 0.0
```
**What:** If the sequence is empty (or was all gaps), return 0.0
instead of crashing with a division-by-zero error.

```python
    gc_count = bases.count("G") + bases.count("C")
    return (gc_count / len(bases)) * 100.0
```
**`.count("G")`** — counts how many times "G" appears in the string.
**`len(bases)`** — total number of characters.
**The math:** (G+C) / total * 100 = GC percentage.

---

### Function 3: sequence_length(sequence)

```python
    return len(sequence.replace("-", ""))
```
**What:** Remove gaps, then count characters. That's it. One line.
**Why a separate function?** Without it, you'd write
`len(seq.replace("-", ""))` everywhere. With the function, you write
`sequence_length(seq)` — clearer and less error-prone.

---

### Function 4: summarize_fasta(filepath)

```python
    lengths = []
    gc_numerator = 0
    total_bases = 0
```
**What:** Initialize collectors.
- `lengths` — will store the length of every sequence.
- `gc_numerator` — running total of G+C bases across ALL sequences.
- `total_bases` — running total of ALL bases across ALL sequences.

```python
    for _header, seq in parse_fasta(filepath):
```
**`_header`** — the underscore prefix is a Python convention meaning
"I need to unpack this value but I won't use it." We only need the
sequence, not the header, for statistics.
**Tuple unpacking:** `parse_fasta` yields `(header, sequence)` tuples.
This line splits each tuple into two separate variables.

```python
        bases = seq.upper().replace("-", "")
        seq_len = len(bases)
        lengths.append(seq_len)
        total_bases += seq_len
        gc_numerator += bases.count("G") + bases.count("C")
```
**`+=`** — add to existing value. `total_bases += 5` is the same as
`total_bases = total_bases + 5`.

```python
    if not lengths:
        raise ValueError(f"No sequences found in {filepath}")
```
**What:** If no sequences were found (empty file), raise an error.
`not lengths` is True when the list is empty.

```python
    overall_gc = (gc_numerator / total_bases) * 100.0 if total_bases else 0.0
```
**Conditional expression:** `A if condition else B`. If total_bases is
not zero, calculate GC%. Otherwise, return 0.0 (avoids division by zero).

```python
    return {
        "file": os.path.basename(filepath),
        "n_seqs": len(lengths),
        "total_bp": total_bases,
        "min_len": min(lengths),
        "mean_len": round(mean(lengths), 1),
        "max_len": max(lengths),
        "gc_percent": round(overall_gc, 2),
    }
```
**What:** Returns a dictionary with all the statistics.
**`min(lengths)`** — smallest value in the list.
**`max(lengths)`** — largest value.
**`round(value, 1)`** — rounds to 1 decimal place. `round(4.667, 1)`
becomes `4.7`.

---

### Function 5: filter_sequences(filepaths, min_len)

```python
    for filepath in filepaths:
        for header, seq in parse_fasta(filepath):
            if sequence_length(seq) >= min_len:
                yield (header, seq)
```
**What:** Nested loops. Outer loop goes through each file. Inner loop
goes through each sequence in that file. If the sequence is long enough,
yield it. Otherwise, skip it (it just doesn't yield).
**`>=`** — greater than or equal to.

---

### Function 6: write_fasta(records, output_path, line_width=80)

```python
def write_fasta(records, output_path, line_width=80):
```
**`line_width=80`** — a default parameter. If you call
`write_fasta(records, "out.fasta")`, line_width is 80. You can override
it: `write_fasta(records, "out.fasta", line_width=60)`.

```python
    with open(output_path, "w") as out_fh:
```
**`"w"`** — write mode (creates or overwrites the file).
**`out_fh`** — short for "output file handle."

```python
        for header, seq in records:
            out_fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                out_fh.write(seq[i:i + line_width] + "\n")
```
**`out_fh.write()`** — writes text to the file (unlike `print`, it
doesn't add a newline, so you add `\n` yourself).
**`range(0, len(seq), line_width)`** — generates numbers 0, 80, 160,
240, etc. This steps through the sequence in chunks of 80.
**`seq[i:i + line_width]`** — slicing. Takes 80 characters starting
at position `i`. This wraps long sequences across multiple lines.

---

### Function 7: format_summary_tsv(summaries)

```python
TSV_COLUMNS = ["file", "n_seqs", "total_bp", "min_len", "mean_len",
               "max_len", "gc_percent"]
```
**What:** A module-level constant (not inside any function). Lists the
column names in order. Defined outside so it can be reused.

```python
    lines = ["\t".join(TSV_COLUMNS)]
```
**`"\t".join(list)`** — joins list items with a tab character between
them. `"\t".join(["a", "b", "c"])` becomes `"a\tb\tc"`.
**`\t`** — the tab character. This is what makes it "tab-separated."

```python
    for summary in summaries:
        row = "\t".join(str(summary[col]) for col in TSV_COLUMNS)
        lines.append(row)
```
**`summary[col]`** — accesses a dictionary value by key.
**`str(...)`** — converts numbers to strings (needed for `.join()`).
**`str(summary[col]) for col in TSV_COLUMNS`** — a generator expression.
It's like a compact for loop inside the `.join()` call.

```python
    return "\n".join(lines) + "\n"
```
**`"\n".join(lines)`** — joins all rows with newlines between them.
**`+ "\n"`** — adds a final newline at the end of the file.

---

### Function 8: build_argument_parser()

```python
    parser = argparse.ArgumentParser(
        description="...",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="...",
    )
```
**`ArgumentParser`** — creates a new command-line parser.
**`description`** — text shown at the top of --help.
**`formatter_class=RawDescriptionHelpFormatter`** — tells argparse to
keep the formatting of the epilog text as-is (don't reflow it).
**`epilog`** — text shown at the bottom of --help (the examples).

```python
    parser.add_argument(
        "fasta_files",
        nargs="+",
        metavar="FASTA",
        help="One or more FASTA files to process.",
    )
```
**`"fasta_files"`** — no dashes means positional argument (required).
**`nargs="+"`** — accepts one or more values. The `+` means "at least
one." Other options: `"*"` (zero or more), `"?"` (zero or one), a
number like `2` (exactly two).
**`metavar="FASTA"`** — what shows up in the help text as the
placeholder name.

```python
    parser.add_argument(
        "--min-len",
        type=int,
        default=None,
        metavar="N",
        help="...",
    )
```
**`"--min-len"`** — double dash means optional flag.
**`type=int`** — auto-converts the string input to an integer.
**`default=None`** — if not provided, the value is None.
Note: argparse converts `--min-len` to `args.min_len` (dashes become
underscores) when you access it.

---

### Function 9: main(argv=None)

```python
def main(argv=None):
    parser = build_argument_parser()
    args = parser.parse_args(argv)
```
**`argv=None`** — when None, argparse reads from `sys.argv` (the actual
command line). When you pass a list (like in tests), it uses that instead.
**`args`** — an object where each argument becomes an attribute.
After parsing `--min-len 90`, you access it as `args.min_len`.

```python
    if args.min_len is not None and args.filtered_out is None:
        parser.error("--min-len requires --filtered-out")
```
**`parser.error()`** — prints an error message and exits the program.
This enforces that --min-len and --filtered-out must be used together.

```python
    summaries = []
    for filepath in args.fasta_files:
        try:
            summary = summarize_fasta(filepath)
            summaries.append(summary)
        except (FileNotFoundError, ValueError) as err:
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
```
**`try/except`** — tries to run the code. If an error of the specified
type occurs, catches it instead of crashing.
**`as err`** — stores the error object in a variable called `err`.
**`sys.exit(1)`** — exits with code 1 (meaning failure).

```python
    tsv_text = format_summary_tsv(summaries)

    if args.out:
        with open(args.out, "w") as tsv_fh:
            tsv_fh.write(tsv_text)
        print(f"Summary written to {args.out}", file=sys.stderr)
    else:
        print(tsv_text, end="")
```
**`if args.out:`** — if the user provided --out, write to that file.
Otherwise print to screen.
**`print(..., end="")`** — normally print adds a newline at the end.
`end=""` suppresses that because the TSV text already ends with `\n`.

```python
if __name__ == "__main__":
    main()
```
**What:** Only run main() when the file is executed directly.
Already explained earlier.

---

## PART 3: THE TEST FILE EXPLAINED

### Modules used in tests

```python
import tempfile
```
**What:** Creates temporary files that are automatically named with
random characters and stored in your system's temp directory.
**Why:** Tests need fake FASTA files to work with. You don't want
to create permanent files that clutter your project.

```python
import unittest
```
**What:** Python's built-in testing framework.

### Key test patterns

```python
class TestParseFasta(unittest.TestCase):
```
**`class`** — defines a group of related tests.
**`unittest.TestCase`** — your class inherits from TestCase, which
gives it all the assertion methods (assertEqual, assertTrue, etc.).

```python
    def test_single_sequence(self):
```
**`test_`** prefix — unittest ONLY runs methods starting with `test_`.
**`self`** — refers to the test instance. Required for all methods
inside a class.

```python
        self.addCleanup(os.remove, path)
```
**What:** Registers a cleanup action. After this test finishes (pass
or fail), Python will call `os.remove(path)` to delete the temp file.
This prevents leftover files.

### Assertion methods used

| Method | What it checks | Example |
|--------|---------------|---------|
| `assertEqual(a, b)` | a equals b exactly | `assertEqual(5, 5)` |
| `assertAlmostEqual(a, b)` | a equals b within rounding | `assertAlmostEqual(50.001, 50.0)` |
| `assertRaises(Error)` | the code raises this error | `assertRaises(ValueError)` |
| `assertIn(a, b)` | a is found inside b | `assertIn("hello", "hello world")` |
| `assertGreater(a, b)` | a > b | `assertGreater(10, 5)` |
| `assertTrue(x)` | x is True | `assertTrue(all(...))` |

### The helper function

```python
def _write_temp_fasta(content):
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    )
    tmp.write(content)
    tmp.close()
    return tmp.name
```
**`_write_temp_fasta`** — the `_` prefix means "private helper, not a
test." It creates a temporary file, writes content to it, and returns
the file path. `delete=False` means don't auto-delete when closed (the
test's addCleanup handles deletion instead).

### The path setup

```python
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
import faxtools
```
**What this does step by step:**
1. `__file__` — the path of THIS test file
   (`/Users/.../biol7180-faxtools/tests/test_faxtools.py`)
2. `os.path.abspath(__file__)` — makes it absolute
3. `os.path.dirname(...)` — goes up one directory (to `tests/`)
4. `os.path.dirname(...)` again — goes up again (to `biol7180-faxtools/`)
5. `sys.path.insert(0, REPO_ROOT)` — tells Python "look in this directory
   when importing modules"
6. `import faxtools` — now Python can find faxtools.py

**Why:** When you run tests from different directories, Python might not
know where faxtools.py is. This ensures it always finds it.

---

## PART 4: OTHER PYTHON CONCEPTS IN YOUR CODE

### Tuples
```python
yield (header, sequence)
```
A tuple is like a list but you can't change it after creation. The
parentheses are optional — `yield header, sequence` works too. Tuples
are used when you want to return multiple values from a function.

### List comprehension vs generator expression
```python
# Generator expression (in format_summary_tsv)
"\t".join(str(summary[col]) for col in TSV_COLUMNS)

# List comprehension (in tests)
lengths = [faxtools.sequence_length(seq) for _, seq in kept]
```
Both are compact loops. The difference: list comprehensions use `[]`
and build a full list in memory. Generator expressions use `()` and
produce items one at a time (more memory efficient).

### Dictionary access
```python
summary["n_seqs"]       # Access value by key
summary[col]            # col is a variable containing the key name
```

### String formatting with f-strings
```python
f"FASTA file not found: {filepath}"
f"Wrote {len(filtered_records)} sequence(s)"
```
f-strings (formatted string literals) let you embed expressions inside
`{}` brackets. Python evaluates the expression and inserts the result.

### The `all()` function (in tests)
```python
self.assertTrue(all(l >= 8 for l in lengths))
```
`all()` returns True only if EVERY item in the iterable is True.
It's a clean way to check that every length meets the threshold.
