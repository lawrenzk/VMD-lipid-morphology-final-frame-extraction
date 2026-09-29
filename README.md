# VMD-lipid-morphology-final-frame-extraction
It scans your existing directory structure:
lipid_water_assembly/
└── systems/
    ├── bcer/
    │   ├── rep01/
    │   │   ├── production.gro
    │   │   └── production.xtc
    │   ├── rep02/
    │   └── rep03/
    ├── dotap/
    └── ...
    For every replica, it extracts the last frame and generates two molecular structures:
Output	Purpose
.pdb	Lipid-only structure for opening in VMD
.gro	Lipid-only coordinates in GROMACS format
.png	Rendered morphology snapshots
.json	Last frame time, lipid count, PBC information
The PNG contains four different views of the assembly.
3D
Isometric assembly view

XY
Top projection

XZ
Side projection

YZ
Second side projection

These are particularly useful because a bilayer and a vesicle can look similar from certain viewing directions.
## How to run it

First, install the dependencies in your Python environment:
python lipid_morphology_lastframe/extract_lastframe_morphologies.py \
    --root lipid_water_assembly/systems \
    --out lipid_water_assembly/morphology_last_frames \
    --wrapped-preview
    
## Expected output structure
The script automatically organizes the generated files:

lipid_water_assembly/
├── systems/
│
└── morphology_last_frames/
    ├── summary.csv
    │
    ├── bcer/
    │   ├── replica_comparison.png
    │   │
    │   ├── rep01/
    │   │   ├── bcer_rep01_lipid_assembled.pdb
    │   │   ├── bcer_rep01_lipid_wrapped.pdb
    │   │   ├── bcer_rep01_lipid_wrapped.gro
    │   │   ├── bcer_rep01_morphology.png
    │   │   ├── bcer_rep01_wrapped.png
    │   │   ├── metadata.json
    │   │   └── view_in_vmd.tcl
    │   │
    │   ├── rep02/
    │   └── rep03/
    │
    ├── dotap/
    └── ...
