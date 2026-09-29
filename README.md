# Lipid-water final-frame morphology extractor

This script scans `systems/<lipid>/rep*/production.gro` and `production.xtc` and uses the **last saved XTC frame**. It outputs only the lipid beads, never waters or ions. It does not require VMD or `gmx_mpi` for extraction.

## Install (in a user-controlled modern Python/Conda environment)

```bash
python -m pip install 'MDAnalysis>=2.5' numpy scipy matplotlib
python -c 'import MDAnalysis, scipy, matplotlib; print("dependencies OK")'
```

Install into a conda env or user environment, not into the compute node's system Python. On Mike, perform the batch visualization in an interactive compute allocation or submit a small Python-only Slurm job, not a large processing job on the login node. The `gromacs/2021.3` module isn't required to run this Python extractor.

## Run from lipid_screening_project

```bash
python extract_lastframe_morphologies.py \
  --root lipid_water_assembly/systems \
  --out lipid_water_assembly/morphology_last_frames \
  --wrapped-preview
```

If lipid or solvent residue names differ, use e.g. `--lipid-resnames BCER` for a single-lipid folder and `--exclude-resnames W,PW,NA,CL,ION` for batch exclusions. Always verify `metadata.json` lists the expected lipid residue and bead count. Run `--force` to regenerate outputs.

Optional: `--contact-cutoff-nm 1.2` controls contact-based inter-molecule reconstruction. `--tail-regex '^(C[0-9]|D[0-9]|T[0-9])'` only influences **illustrative bead colors**. It does not identify chemistry or classify morphology.

## What you get

```
lipid_water_assembly/morphology_last_frames/
  summary.csv
  bcer/
    replica_comparison.png
    rep01/
      bcer_rep01_lipid_wrapped.pdb    # original box positions; lipid only
      bcer_rep01_lipid_wrapped.gro    # original box positions; lipid only
      bcer_rep01_lipid_assembled.pdb  # view-oriented PBC reconstruction
      bcer_rep01_morphology.png       # 3D and XY/XZ/YZ views
      bcer_rep01_wrapped.png          # only if --wrapped-preview
      metadata.json                  # actual last saved time and selection
      view_in_vmd.tcl                 # open the assembled PDB in VMD later
    rep02/...
```

Open the assembled PDB in VMD (`vmd -e view_in_vmd.tcl`) from its directory. PDB is a *single static frame*, not an XTC movie; the wrapped GRO/PDB preserve original atom positions. A static PNG is for rapid screening; do not infer morphology solely from one projection or from one replica.

### PBC caveats

For typical orthorhombic systems, the extractor makes each lipid residue whole, connects lipids whose beads are within the contact radius under periodic boundaries, then attempts to reassemble each connected component. A periodic planar bilayer may percolate across the box, so this reconstruction can be **non-unique**; `metadata.json` records conflicting periodic graph edges. A tilted/triclinic box is only reconstructed *within individual lipids*, and the metadata flags that limitation. Always keep the original wrapped coordinates; for any uncertain morphology, inspect more frames with explicit topology/bond-aware reconstruction. This script **does not** automatically label a structure as a vesicle or bilayer.

The script assumes the GRO topology contains one lipid per residue and that atom ordering matches XTC. It reports and skips files that are missing or fail to open, and uses actual XTC last frame time (not an assumed 750 ns).
