# ms2Topo integration validation

Validated on 2026-11-21 against OrbiFragsNets base commit
`44063e4b7a7fb2b0af4e28a04d9f6074be0dcc83`, using ms2Topo exports from
`58a4aa48cd0f387eb1e2bc27912562581ca39dd2`.

## Checks

- `python -m pytest -q`: **27 passed**.
- Tests include actual formula/network annotation of exact-mass synthetic
  fragments, conversion of the unchanged ms2Topo CSV, input validation,
  intensity filtering/renormalization, feature-ID matching, no-result cases,
  explicit batch errors, search limits, and use from another working directory.
- Both package imports and the original flat `Functions` imports work for
  annotation without pyOpenMS, IPython, or tabulate installed.
- Vectorized exact masses match the original scalar calculation within
  1e-12 Da for 100 deterministic test formulas. Candidate mass calculation
  reads MassVec.csv once per batch instead of once per formula.
- The example notebook's Python cells parse successfully. The API/CLI calls
  it demonstrates were exercised separately.
- The direct batch CLI completed on real ms2Topo example features 0 and 1
  with default options and the repository's unchanged elemental bounds.

| Feature | Input peaks | Above precursor | Selected usable peaks | Annotated peaks | Status |
|---|---:|---:|---:|---:|---|
| 0 | 32 | 4 | 10 | 8 | annotated |
| 1 | 8 | 2 | 1 | 0 | too_few_fragments |

Selected usable peaks also exclude nonpositive median intensity. These
counts therefore do not equal input peaks minus above-precursor peaks.
See `examples/ms2topo/results/annotations.csv` and `diagnostics.csv` for
the example outputs. Formulas are algorithmic annotations, not validated
compound identifications.

Test environment: Python 3.14.8, NumPy 2.5.5, pandas 3.1.1, SciPy 1.18.1,
pytest 9.2.0. Older dependency versions listed in requirements-ms2topo.txt
were not separately tested.

## Scope and limitations

The adapter starts after consensus formation. It does not rerun raw profile
processing or change the original formula-selection/network-scoring rules.
Consensus IQR is used as a configurable spread/tolerance proxy; this is not
the same uncertainty model as the original raw-profile confidence intervals.
The original positive-ion mass convention and valence constraints remain.
The full raw mzML workflow was not rerun because pyOpenMS was not installed.

Large peak/formula-network searches stop at the configured max_networks bound
rather than enumerating beyond it. Runtime still depends on precursor
candidate count and elemental bounds. The limit bounds network enumeration,
not precursor formula enumeration or total wall-clock time.

All changes remain local on branch `feat/ms2topo-input`.

## Synchronization with updated main

Rebased onto upstream main commit `2471529f1132ede58d4fc5a1811c95f1624f2365`.
The upstream batch runner, adapter, documentation, datasets, and updated
MassVec/MaxAtomicSubscripts parameter files are preserved.

After synchronization, all **28 tests passed**, including an end-to-end test
of the existing upstream compact-formula batch runner. Checks used a temporary
copy of the same code, current parameters, and small fixtures, excluding the
large dataset and wheel archive directories. The example annotations above
were generated before this parameter update and remain historical examples.
