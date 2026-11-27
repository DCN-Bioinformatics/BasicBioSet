# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Python 3.10+, standard library only — there is no requirements file, build step, or linter config.

```bash
python3 -m unittest discover tests                                   # all tests
python3 -m unittest tests.test_sequence                              # one module
python3 -m unittest tests.test_web.TestLoadText.test_fasta           # one test

python3 main.py                      # interactive menu on data/example.fasta
python3 main.py file.fasta           # full report, no menu (--min-orf N to filter)
python3 main.py file.fasta --menu    # menu with a custom file
python3 main.py --web [--port 8080]  # serve the website at http://127.0.0.1:8000/web/
```

The website must be served over HTTP (via `--web` or GitHub Pages); opening `web/index.html` from disk fails because the browser blocks fetching the `.py` files.

## Architecture

One analysis library with two front ends that must give identical results:

- **`bioseq/`** — the pure-Python core (FASTA parsing, sequence ops, motifs/restriction enzymes/k-mers, Tm and molecular weight, sliding-window composition and CpG islands). FASTQ reading, quality stats and trimming live in `fasta.py` (deliberately not a new module: see the caching note below). `validate()` in `sequence.py` is the single entry point for input: every public function calls it, and it uppercases and reads RNA U as T. Everything public is re-exported from `bioseq/__init__.py`; both front ends import from the package root.
- **`main.py`** — terminal UI. Report mode, an interactive menu (`mode_*` functions, one per menu entry), and `--web`, which starts a static `http.server` rooted at the repo so `web/` can reach `../bioseq/` and `../data/`.
- **`web/`** — browser UI that runs the *real* `bioseq` package inside Pyodide (loaded from jsDelivr, pinned in `web/app.js`). No server-side analysis; sequences never leave the browser.
  - `web/app.js` fetches the files listed in `BIOSEQ_FILES` plus `bridge.py`, writes them into Pyodide's virtual FS under `/app`, and imports `bridge`. **Adding a new module to `bioseq/` requires adding it to `BIOSEQ_FILES`**, or the web app will fail to import it.
  - Browsers cache `app.js` (10 min on GitHub Pages) but fetch the `.py` files with `cache: "no-cache"`, so right after a deploy a stale `BIOSEQ_FILES` can meet a new `__init__.py` that imports a module missing from the list, and the engine fails to start. Prefer extending existing modules over adding new ones.
  - `web/bridge.py` is the JSON boundary: every public function takes strings and returns a JSON string `{"ok": true, ...}` or `{"ok": false, "error": "..."}`, so JS never touches Python objects. `load_text()` parses input into module-level `_records`; `run(tool, seq_id, params_json)` dispatches through the `TOOLS` dict; `report` and `quality` are special-cased because they cover all loaded records (`DATASET_TABS` in `app.js`). FASTQ input keeps per-read Phred scores in `_qualities` alongside `_records`. A new web analysis = a `_tool(seq, params) -> dict` function registered in `TOOLS`, plus rendering code in `app.js`.
  - Positions: `bioseq.find_orfs` returns 0-based `start`; the bridge and terminal output convert to 1-based for display.
- **`tests/`** — `test_sequence.py` covers the library; `test_web.py` imports `web/bridge.py` in regular CPython (it runs unchanged outside Pyodide) and asserts on the decoded JSON.

## Constraints

- Keep `bioseq` and `bridge.py` dependency-free: they run in Pyodide, where only files fetched by `app.js` and the standard library are available.
- The empty `.nojekyll` at the repo root is required for GitHub Pages, which otherwise skips `bioseq/__init__.py`.
- `.github/workflows/`: `tests.yml` runs the unit tests on Python 3.10 and 3.12 for every PR and push to `main`; the Claude workflows review every PR and respond to `@claude` mentions in issues and PRs.
