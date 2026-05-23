# DNA Sequence Assembly

**Advanced Programming in Python — Project**

| | |
|---|---|
| **Author** | Ihame Gilbert |

Program that assembles DNA from overlapping segments using a **de Bruijn graph** and an **Euler path**.

---

## Requirements

- Python 3.10+
- See `requirements.txt` for package versions

### Install dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install pandas networkx matplotlib pytest
```

---

## Project structure

```
Project/
├── project.py          # Main implementation and entry point
├── project_test.py     # Unit tests (pytest)
├── requirements.txt    # Python dependencies
├── README.md           # This file
├── csv-files/          # Input CSV files (local development)
│   ├── DNA_1_5.csv
│   ├── DNA_2_5.csv
│   ├── DNA_3_3.csv
│   ├── DNA_4_5.csv
│   └── DNA_5_5.csv
├── DNA_x.png           # Output: graph plot (generated)
└── DNA_x.txt           # Output: DNA sequence (generated)
```

> **Submission note:** The assignment requires submitting files in a single folder without subdirectories. For local development, CSV files may remain in `csv-files/`.

---

## Run the program

Open a terminal in the project folder:

```bash
cd C:/Users/gihame/Music/Project
```

Run with a CSV file named `DNA_x_k.csv` (`x` = dataset id, `k` = k-mer length):

```bash
python project.py csv-files/DNA_1_5.csv
```

Other examples:

```bash
python project.py csv-files/DNA_2_5.csv
python project.py csv-files/DNA_3_3.csv
python project.py csv-files/DNA_4_5.csv
python project.py csv-files/DNA_5_5.csv
```

### Expected output

1. **Console** — Euler path printed by `construct_dna_sequence` (e.g. `CTGA -> TGAA -> GAAT -> ...`)
2. **`DNA_x.png`** — de Bruijn graph visualization
3. **`DNA_x.txt`** — reconstructed DNA sequence, or `DNA sequence can not be constructed.`

Example: `DNA_1_5.csv` creates `DNA_1.png` and `DNA_1.txt`.

### View results

**DNA sequence:**

```bash
cat DNA_1.txt
```

**Graph image:**

- Open `DNA_x.png` in VS Code / Cursor or double-click in File Explorer
- Windows (Git Bash): `start DNA_1.png`

**Optional pop-up graph:** Add `plt.show()` before `plt.close()` in `plot_graph` to display the graph in a window when running.

---

## Pipeline

| Task | Function | Description |
|------|----------|-------------|
| 1 | `read_csv` | Load CSV into pandas DataFrame |
| 2 | `clean_data` | Remove invalid/duplicate segments |
| 3 | `generate_sequences` | Build JSON of segment sequences |
| 4 | `construct_graph` | Build de Bruijn graph (NetworkX MultiDiGraph) |
| 5 | `plot_graph` | Save graph as `DNA_x.png` |
| 6 | `is_valid_graph` | Check Euler path conditions |
| 7 | `construct_dna_sequence` | Build DNA string and print Euler path |
| 8 | `save_output` | Write result to `DNA_x.txt` |

`main()` parses `DNA_x_k.csv` from the command line and runs the full pipeline.

---

## Run tests

From the project folder:

```bash
pytest project_test.py -v
```

Individual tests:

```bash
pytest project_test.py::test_clean_data -v
pytest project_test.py::test_generate_sequences -v
pytest project_test.py::test_construct_graph -v
pytest project_test.py::test_is_valid_graph -v
pytest project_test.py::test_construct_dna_sequence -v
```

Tasks **1, 5, and 8** are verified by running `main` and checking output files (no required pytest).

---

## Quick manual checks

**Read CSV:**

```bash
python -c "from project import read_csv; df=read_csv('csv-files/DNA_1_5.csv'); print(df.head()); print('rows:', len(df))"
```

**Generate JSON:**

```bash
python -c "from project import read_csv, clean_data, generate_sequences; df=clean_data(read_csv('csv-files/DNA_1_5.csv')); print(generate_sequences(df))"
```

**Small graph example:**

```bash
python -c "import json; from project import construct_graph; G=construct_graph(json.dumps({'1':'ATTACTC'}),5); print(sorted(G.edges()))"
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Invalid filename format` | Filename must match `DNA_x_k.csv`; `main` uses `os.path.basename` for paths like `csv-files/...` |
| `NameError: os` | Add `import os` in `project.py` |
| `NameError: save_output` | Define `save_output` in `project.py` |
| `NameError: json` | Add `import json` in `project.py` |
| `FileNotFoundError` | Run from the project folder; use `csv-files/DNA_1_5.csv` |
| No `DNA_x.png` | Check terminal for errors before `plot_graph` |

---

## Bonus tasks (optional)

- **Task 9:** `construct_all_dna` → `DNA_x_All.txt`
- **Task 10:** `align_segments` → `DNA_x_alignment_<dna>.txt`
- **Task 11:** `generate_all_alignments`

---

## Project report

- **`REPORT.md`** — Draft report (convert to PDF for submission)
- **`generate_report_outputs.py`** — Prints intermediate results for the example run section

```bash
python generate_report_outputs.py
```

## References

Course project: *DNA Sequence Assembly* (de Bruijn graph, Euler path, Sanger CSV format).
