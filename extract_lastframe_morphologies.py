#!/usr/bin/env python3
"""Extract and visualize the final lipid-only frame from GROMACS GRO+XTC replicas.

Requires numpy, scipy, matplotlib, MDAnalysis. Read README.md for details.
The PBC reconstruction is a visualization heuristic, not a substitute for
bond-based whole-molecule reconstruction or a morphology classifier.
"""
import argparse
import csv
import json
import re
import sys
import warnings
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.mdamath import triclinic_vectors
    from scipy.spatial import cKDTree
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError as exc:
    sys.exit('Missing Python dependency: {}\nInstall: python -m pip install MDAnalysis scipy matplotlib numpy'.format(exc))

DEFAULT_EXCLUDE = 'W,PW,WF,SOL,HOH,WAT,TIP3,TIP3P,NA,NA+,CL,CL-,K,K+,CA,CA2+,MG,MG2+,ION,SOD,CLA'
DEFAULT_TAIL = r'^(?:C\d|D\d|T\d)'


def names_csv(arg):
    return {x.strip().upper() for x in arg.split(',') if x.strip()}


def unitcell(ts):
    box = np.asarray(ts.dimensions, dtype=float)
    if box.shape[0] != 6 or min(box[:3]) <= 0:
        raise ValueError('Frame lacks a valid periodic box')
    cell = np.asarray(triclinic_vectors(box), dtype=float)
    return box, cell


def residue_groups(atomgroup):
    """Indices in the selected atom group's coordinate array for each molecule.

    This assumes one lipid per GRO residue, as in the provided Packmol workflow.
    """
    chosen = atomgroup.indices
    groups = []
    for residue in atomgroup.residues:
        idx = np.searchsorted(chosen, residue.atoms.indices)
        if not np.array_equal(chosen[idx], residue.atoms.indices):
            raise ValueError('Incomplete lipid residue encountered; check residue selection')
        groups.append(idx)
    return groups


def make_each_lipid_whole(xyz, groups, cell):
    """Minimum-image reconstruct relative to first bead of each lipid (Å)."""
    invcell = np.linalg.inv(cell)
    whole = xyz.copy()
    for ids in groups:
        anchor = xyz[ids[0]]
        dr_fraction = (xyz[ids] - anchor) @ invcell
        dr_fraction -= np.rint(dr_fraction)
        whole[ids] = anchor + dr_fraction @ cell
    return whole


def assemble_orthorhombic(whole, groups, lengths, contact_angstrom):
    """Build a PBC lipid-contact graph; unwrap each contact-connected component.

    A periodic bilayer can percolate, producing contradictory graph cycles. Such
    cases are flagged; always consult the saved wrapped structure alongside it.
    """
    n = len(groups)
    if n == 0:
        raise ValueError('No lipid residues selected')
    atom_to_res = np.empty(len(whole), dtype=np.int32)
    centers = np.empty((n, 3), dtype=float)
    for i, ids in enumerate(groups):
        atom_to_res[ids] = i
        centers[i] = whole[ids].mean(axis=0)

    adj = [set() for _ in range(n)]
    tree = cKDTree(np.mod(whole, lengths), boxsize=lengths)
    pairs = tree.query_pairs(contact_angstrom, output_type='ndarray')
    for a, b in pairs:
        i, j = int(atom_to_res[a]), int(atom_to_res[b])
        if i != j:
            adj[i].add(j)
            adj[j].add(i)

    shifts = np.zeros((n, 3), dtype=float)
    visited = set()
    components = []
    conflict_edges = set()
    for root in range(n):
        if root in visited:
            continue
        visited.add(root)
        queue = deque([root])
        component = []
        while queue:
            i = queue.popleft()
            component.append(i)
            for j in sorted(adj[i]):
                delta = centers[j] - centers[i]
                min_delta = delta - lengths * np.rint(delta / lengths)
                expected = shifts[i] + min_delta - delta
                if j not in visited:
                    shifts[j] = expected
                    visited.add(j)
                    queue.append(j)
                elif np.max(np.abs(shifts[j] - expected)) > 1.0:
                    conflict_edges.add(tuple(sorted((i, j))))
        components.append(component)

    # Keep independent assemblies near the largest assembly, without merging them.
    largest = max(components, key=len)
    primary_center = (centers[largest] + shifts[largest]).mean(axis=0)
    for component in components:
        if component == largest:
            continue
        centroid = (centers[component] + shifts[component]).mean(axis=0)
        delta = centroid - primary_center
        component_shift = -lengths * np.rint(delta / lengths)
        shifts[component] += component_shift

    reconstructed = whole.copy()
    for i, ids in enumerate(groups):
        reconstructed[ids] += shifts[i]

    # Center the largest assembly for consistent presentation.
    largest_atoms = np.concatenate([groups[i] for i in largest])
    move = lengths / 2 - reconstructed[largest_atoms].mean(axis=0)
    reconstructed += move
    return reconstructed, sorted((len(c) for c in components), reverse=True), len(conflict_edges)


def style_axes_3d(ax, xyz_nm, title, span):
    middle = (xyz_nm.max(axis=0) + xyz_nm.min(axis=0)) / 2
    half = span / 2
    ax.set_xlim(middle[0] - half, middle[0] + half)
    ax.set_ylim(middle[1] - half, middle[1] + half)
    ax.set_zlim(middle[2] - half, middle[2] + half)
    ax.set_xlabel('X (nm)', fontsize=8)
    ax.set_ylabel('Y (nm)', fontsize=8)
    ax.set_zlabel('Z (nm)', fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_box_aspect((1, 1, 1))
    ax.tick_params(labelsize=7)
    ax.view_init(elev=22, azim=48)


def plot_snapshot(xyz_angstrom, names, title, dest, tail_pattern, dpi):
    xyz_nm = xyz_angstrom / 10.0
    if len(xyz_nm) == 0:
        raise ValueError('No atoms for rendering')
    is_tail = np.array([bool(tail_pattern.search(str(n))) for n in names])
    colors = np.where(is_tail, '#168B88', '#D65B35')
    spread = np.ptp(xyz_nm, axis=0)
    span = max(float(spread.max()) * 1.12, 1.0)
    point_size = max(1.0, min(12.0, 18000.0 / len(xyz_nm)))

    fig = plt.figure(figsize=(11, 9), facecolor='white')
    fig.subplots_adjust(left=.08, right=.96, bottom=.10, top=.90, hspace=.36, wspace=.26)
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.scatter(xyz_nm[:, 0], xyz_nm[:, 1], xyz_nm[:, 2], c=colors,
               s=point_size, alpha=0.72, depthshade=True, linewidths=0, rasterized=True)
    style_axes_3d(ax, xyz_nm, 'Isometric', span)
    for panel, (ia, ib, label) in enumerate(((0, 1, 'XY / top'), (0, 2, 'XZ / side'), (1, 2, 'YZ / side')), start=2):
        ax2 = fig.add_subplot(2, 2, panel)
        ax2.scatter(xyz_nm[is_tail, ia], xyz_nm[is_tail, ib], s=point_size,
                    c='#168B88', alpha=.48, linewidths=0, rasterized=True, label='Tail-name beads')
        ax2.scatter(xyz_nm[~is_tail, ia], xyz_nm[~is_tail, ib], s=point_size,
                    c='#D65B35', alpha=.72, linewidths=0, rasterized=True, label='Other beads')
        middle = (xyz_nm[:, [ia, ib]].max(axis=0) + xyz_nm[:, [ia, ib]].min(axis=0)) / 2
        ax2.set_xlim(middle[0]-span/2, middle[0]+span/2)
        ax2.set_ylim(middle[1]-span/2, middle[1]+span/2)
        ax2.set_aspect('equal', adjustable='box')
        ax2.set_xlabel(('X', 'Y', 'Z')[ia] + ' (nm)')
        ax2.set_ylabel(('X', 'Y', 'Z')[ib] + ' (nm)')
        ax2.set_title(label)
        ax2.grid(alpha=.14)
        ax2.tick_params(labelsize=8)
    fig.suptitle(title, fontsize=13, weight='bold', y=.975)
    fig.text(.5, .028, 'Illustrative atom-name coloring (not chemically validated): orange = other, teal = tail-name beads.',
             ha='center', fontsize=8)
    fig.savefig(dest, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def write_vmd_script(outdir, pdb_name):
    tcl = '''# Run locally from this folder: vmd -e view_in_vmd.tcl
mol new {%s} type pdb waitfor all
mol delrep 0 top
mol selection all
mol representation VDW 0.9 12
mol color Name
mol material Opaque
mol addrep top
display projection Orthographic
color Display Background white
display resetview
''' % pdb_name
    (outdir / 'view_in_vmd.tcl').write_text(tcl)


def analyze_replica(rep, output, args):
    gro = rep / args.gro_name
    xtc = rep / args.xtc_name
    if not gro.is_file() or not xtc.is_file() or xtc.stat().st_size == 0:
        missing = ', '.join(str(p.name) for p in (gro, xtc) if not p.is_file() or p.stat().st_size == 0)
        raise FileNotFoundError('Missing/empty input: {}'.format(missing))

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='.*element information.*', category=UserWarning)
        u = mda.Universe(str(gro), str(xtc))
        ts = u.trajectory[-1]  # last *saved* frame, NOT assumed 750 ns
    frame_num, time_ps = int(ts.frame), float(ts.time)
    box, cell = unitcell(ts)
    requested = names_csv(args.lipid_resnames) if args.lipid_resnames else None
    excluded = names_csv(args.exclude_resnames)
    all_names = np.array([str(n).upper() for n in u.atoms.resnames])
    mask = np.isin(all_names, list(requested)) if requested else ~np.isin(all_names, list(excluded))
    selected = u.atoms[np.flatnonzero(mask)]
    if len(selected) == 0:
        raise ValueError('No lipid atoms identified. Use --lipid-resnames NAME.')
    groups = residue_groups(selected)
    xyz = selected.positions.astype(np.float64).copy()  # MDAnalysis uses Å
    if not np.isfinite(xyz).all():
        raise ValueError('Non-finite final-frame coordinates')
    outdir = output
    outdir.mkdir(parents=True, exist_ok=True)
    base = '{}_{}'.format(rep.parent.name, rep.name)
    wrapped_pdb = outdir / (base + '_lipid_wrapped.pdb')
    wrapped_gro = outdir / (base + '_lipid_wrapped.gro')
    assembled_pdb = outdir / (base + '_lipid_assembled.pdb')
    png = outdir / (base + '_morphology.png')
    if not args.force and all(x.is_file() for x in (assembled_pdb, wrapped_gro, png, outdir / 'metadata.json')):
        with (outdir / 'metadata.json').open() as fh:
            return json.load(fh), 'skipped_existing'

    # Save original simulation coordinates unchanged before making a view-oriented copy.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        selected.write(str(wrapped_pdb))
        selected.write(str(wrapped_gro))
    whole = make_each_lipid_whole(xyz, groups, cell)
    if np.max(np.abs(box[3:] - 90)) > 0.05:
        # Do not pretend that an orthorhombic periodic KDTree handles tilted cells.
        reconstructed = whole.copy()
        component_sizes, conflicts = [len(groups)], -1
        pbc_note = 'Triclinic box: molecules made whole; inter-molecular assembly not reconstructed.'
    else:
        reconstructed, component_sizes, conflicts = assemble_orthorhombic(
            whole, groups, box[:3], args.contact_cutoff_nm * 10.0)
        pbc_note = ('PBC cycle conflicts detected: assembled view is approximate; inspect wrapped PDB.'
                    if conflicts else 'Contact-graph PBC reconstruction for visualization only.')

    selected.positions = reconstructed.astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        selected.write(str(assembled_pdb))
    tail_pattern = re.compile(args.tail_regex, re.IGNORECASE)
    plot_snapshot(reconstructed, selected.names, '{} · {} · last saved frame {:.3f} ns'.format(
        rep.parent.name, rep.name, time_ps/1000), png, tail_pattern, args.dpi)
    if args.wrapped_preview:
        plot_snapshot(xyz, selected.names, '{} · {} · wrapped / {:.3f} ns'.format(
            rep.parent.name, rep.name, time_ps/1000), outdir / (base+'_wrapped.png'), tail_pattern, args.dpi)
    write_vmd_script(outdir, assembled_pdb.name)
    meta = {
        'lipid_directory': rep.parent.name, 'replica': rep.name,
        'frame_index_zero_based': frame_num, 'last_saved_time_ps': time_ps,
        'last_saved_time_ns': time_ps / 1000.0,
        'lipid_resnames': dict(Counter(map(str, selected.residues.resnames))),
        'lipid_residues': len(groups), 'lipid_beads': len(selected),
        'selected_atom_resnames': sorted(set(map(str, selected.resnames))),
        'box_dimensions_angstrom': [float(x) for x in box],
        'component_sizes_lipids': component_sizes, 'pbc_cycle_conflicts': conflicts,
        'pbc_note': pbc_note, 'tail_regex': args.tail_regex,
        'wrapped_pdb': wrapped_pdb.name, 'wrapped_gro': wrapped_gro.name,
        'assembled_pdb': assembled_pdb.name, 'snapshot_png': png.name,
        'interpretation': 'Visualization only; no automated morphology classification.'
    }
    with (outdir / 'metadata.json').open('w') as fh:
        json.dump(meta, fh, indent=2)
    return meta, 'ok'


def contact_sheets(outroot, dpi):
    for lipid_dir in sorted(outroot.iterdir()):
        if not lipid_dir.is_dir():
            continue
        images = []
        for repdir in sorted(lipid_dir.glob('rep*')):
            imgs = sorted(repdir.glob('*_morphology.png'))
            if imgs:
                images.append((repdir.name, imgs[0]))
        if len(images) < 2:
            continue
        fig, axs = plt.subplots(1, len(images), figsize=(5*len(images), 4.5), squeeze=False)
        for ax, (label, p) in zip(axs[0], images):
            ax.imshow(plt.imread(p))
            ax.set_title(label)
            ax.axis('off')
        fig.suptitle(lipid_dir.name + ' – last frame by replica')
        fig.savefig(lipid_dir / 'replica_comparison.png', dpi=min(160,dpi), bbox_inches='tight')
        plt.close(fig)


def get_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--root', type=Path, default=Path('lipid_water_assembly/systems'), help='systems/ directory')
    p.add_argument('--out', type=Path, default=None, help='central result directory')
    p.add_argument('--gro-name', default='production.gro')
    p.add_argument('--xtc-name', default='production.xtc')
    p.add_argument('--lipid-resnames', default='', help='Optional comma-separated allow-list, e.g. BCER,DPPC')
    p.add_argument('--exclude-resnames', default=DEFAULT_EXCLUDE, help='Comma-separated solvents/ions to exclude')
    p.add_argument('--contact-cutoff-nm', type=float, default=1.2, help='Bead-bead contact graph radius for view-oriented PBC reconstruction')
    p.add_argument('--tail-regex', default=DEFAULT_TAIL, help='Atom-name regex used only for preview color coding')
    p.add_argument('--dpi', type=int, default=170)
    p.add_argument('--wrapped-preview', action='store_true', help='Also plot a wrapped-coordinate PNG')
    p.add_argument('--force', action='store_true', help='Overwrite previously created results')
    p.add_argument('--no-comparisons', action='store_true', help='Skip replica-comparison contact sheets')
    return p.parse_args(argv)


def main(argv=None):
    args = get_args(argv)
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit('No such systems directory: {} (use --root)'.format(root))
    outroot = (args.out or root.parent / 'morphology_last_frames').expanduser().resolve()
    if outroot == root or root in outroot.parents:
        raise SystemExit('Output must NOT be placed inside systems/; choose another --out directory')
    outroot.mkdir(parents=True, exist_ok=True)
    dirs = sorted({p.parent for patt in (args.gro_name, args.xtc_name)
                   for p in root.rglob(patt) if p.parent != root})
    if not dirs:
        raise SystemExit('No replica folders found containing {} / {}'.format(args.gro_name,args.xtc_name))
    print('Found {} replica directories; output {}'.format(len(dirs),outroot), flush=True)
    results = []
    for k, rep in enumerate(dirs, 1):
        label = str(rep.relative_to(root))
        output = outroot / rep.relative_to(root)
        try:
            meta, status = analyze_replica(rep, output, args)
            print('[{:d}/{:d}] {:<45} {} {:.3f} ns ({} lipid molecules, {} PBC conflicts)'.format(
                k, len(dirs), label, status, meta['last_saved_time_ns'],
                meta['lipid_residues'], meta['pbc_cycle_conflicts']), flush=True)
            results.append({'system': label, 'status': status, 'time_ns': meta['last_saved_time_ns'],
                            'n_lipids': meta['lipid_residues'], 'n_beads': meta['lipid_beads'],
                            'components': ','.join(map(str,meta['component_sizes_lipids'])),
                            'pbc_conflicts': meta['pbc_cycle_conflicts'], 'detail': meta['pbc_note']})
        except Exception as exc:
            print('[{:d}/{:d}] {} ERROR: {}'.format(k,len(dirs),label,exc), file=sys.stderr, flush=True)
            results.append({'system': label, 'status':'error', 'time_ns':'', 'n_lipids':'', 'n_beads':'',
                            'components':'', 'pbc_conflicts':'', 'detail':str(exc)})
    fields = ['system','status','time_ns','n_lipids','n_beads','components','pbc_conflicts','detail']
    with (outroot/'summary.csv').open('w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(results)
    if not args.no_comparisons:
        contact_sheets(outroot, args.dpi)
    ok = sum(r['status'] != 'error' for r in results)
    print('Completed {}/{}; see {}'.format(ok,len(results),outroot/'summary.csv'))
    return 0 if ok else 2


if __name__ == '__main__':
    sys.exit(main())
