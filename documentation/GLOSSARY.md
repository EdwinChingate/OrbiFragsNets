# Glossary

- **m/z**: Mass-to-charge ratio measured by the Orbitrap analyzer; the primary coordinate for peaks and annotations.
- **MS1 / MS2**: MS1 scans capture precursor ion profiles across RT; MS2 scans isolate a precursor and record its fragmentation spectrum.
- **Orbitrap**: Electrostatic ion trap where ion oscillations are converted into high-resolution m/z measurements via Fourier transform.
- **Centroiding**: Conversion of raw profile spectra into discrete peaks (m/z, intensity) to reduce data volume while retaining peak centroids.
- **Profile mode**: Raw acquisition format containing the full time-domain signal before centroiding.
- **Retention time (RT)**: Chromatographic time coordinate describing when a compound elutes; used to align MS1 and MS2 events.
- **ppm (parts per million)**: Relative error unit used for mass differences and confidence intervals.
- **Fragment**: Product ion observed in MS2 after the precursor ion undergoes fragmentation.
- **Annotation**: Assignment of a molecular formula (and sometimes adduct) to a fragment based on exact mass, intensity, and network consistency.
- **Chemical space**: Set of candidate molecular formulas that satisfy the configured mass and element-count constraints.
- **Chemical consistency**: Requirement that elemental balances between fragments and neutral losses are satisfied, modeled as vector equations on element counts.
- **Fragments network**: Graph whose nodes are fragment annotations and whose edges represent chemically consistent relationships or neutral losses.
- **Network grade**: Numerical score summarizing how well-connected a fragment network is; used to rank competing annotation sets.
- **DDA (data-dependent acquisition)**: Acquisition mode in which MS2 scans are triggered dynamically based on MS1 signals.
