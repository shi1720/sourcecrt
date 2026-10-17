# SourceCRT

Exact source-aware certification of bounded integer vectors from possibly corrupted modular images.

**Manuscript:** Shivam Gupta, *Source-Aware Certificates for Bounded Modular Reconstruction* (2026), an unpublished research manuscript. See `paper/main.pdf`. This repository does not imply acceptance, peer review, or a new general CRT decoding algorithm.

## Contract

Inputs are pairwise coprime moduli, received integer vectors, trusted integer intervals, and a finite family of allowed corruption supports. At least one support must contain all erroneous image positions. Source groups can overlap. The method exactly returns inconsistent, unique, or ambiguous under this declaration. It cannot validate whether a physical computation obeyed that declaration.

The generator uses incremental CRT; the separate verifier uses a direct CRT sum and adjacent progression checks. Core package has no external dependencies. Verification checks remain active under Python `-O`.

## Reproduce

Python 3.10+; experiments originally run on Python 3.12.14. Install in a virtual environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[experiments]'
python -m unittest discover -s tests -v
python -O -m unittest discover -s tests -v
python experiments/run.py
python experiments/audit.py
python experiments/make_tables.py
cd paper
tectonic main.tex
```

Tectonic is only required to compile the PDF. Seeds, objects, bounds, prime ordering, and fault grid are fixed in `experiments/run.py`. Experimental dependencies are pinned in `pyproject.toml`; this is not a fully locked operating-system image. Generation scripts use normal Python assertions for research controls, so run them without `-O`.

## Use

```python
from sourcecrt import certificate, verify

moduli = [5, 7, 11]
residues = [[1, 2], [0, 0], [0, 0]]
bounds = [(0, 11), (0, 11)]
supports = [(0,), (1,), (2,)]  # one faulty source, one image each
cert = certificate(moduli, residues, bounds, supports)
assert cert['candidate'] == [0, 0]
assert verify(moduli, residues, bounds, supports, cert)
```

Indices are zero-based in code and one-based in the paper. The two scalar problems in this example are both ambiguous. No denominator reconstruction or polynomial coefficient relations are assumed. This is a small-support-family tool: enumeration grows with the number of admissible source unions.

## Evidence

- 48 norm polynomials in fields defined by `t**d - t - 1`, for degrees 4, 8, 12, 16.
- Exact characteristic-polynomial and resultant ground truth independently agree.
- 2,304 modular images from a separate finite-field Faddeev-LeVerrier implementation.
- 576 fault/placement/budget configurations, with all original transcripts and certificates.
- 4,500 exhaustive small received arrays checked against direct point enumeration.
- 1,728 certificate mutations rejected; 1,200 exact placement checks.
- An explicit out-of-model common-mode counterexample is included.

`results/trials.csv`: all stopping counts and diagnostic timings. `objects.json`: input elements, multiplication matrices, exact norms and bounds. `transcripts.json`: all full received arrays, source owners and declared fault IDs. `certificates.json`: first joint stopping declaration and certificate for each configuration. `audit.json` and `summary.json`: aggregate checks. Timings vary on reruns; integer counts and mathematical outputs are deterministic for the stated environment. The 576 configurations reuse 48 objects and are not independent samples.

Joint exact stopping saves prime computations over exact scalar stopping in 182/576 configurations, with a maximum saving of six. The single-coefficient error control has no saving. These are evidence counts, not cluster runtime speedups. Existing interleaved lattice decoders are discussed and credited in the paper; no implementation-level comparison against those decoders is claimed.

## Prior art

Error-tolerant CRT, rational reconstruction, and interleaved CRT decoding are established. Especially relevant are Goldreich–Ron–Sudan (2000), Böhm–Decker–Fieker–Pfister (2015), Abbott (2017), and Li–Sidorenko–Nielsen (ISIT 2013). Full references and the limited contribution claim are in the manuscript.

## Layout

- `src/sourcecrt/core.py`: progression boxes, generation, universal bound, source layout
- `src/sourcecrt/verify.py`: independent verifier
- `tests/test_exact.py`: direct exhaustive and edge-case controls
- `experiments/`: deterministic reproduction and audit scripts
- `results/`: full evidence, not only aggregate tables
- `figures/`: vector PDF figures generated from the evidence
- `paper/`: editable LaTeX source, generated tables, manuscript PDF

Code: MIT license. Manuscript text and figures: copyright Shivam Gupta, 2026; all rights reserved. No third-party papers or private profile documents are redistributed.
