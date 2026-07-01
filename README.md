# MolScreener

**Compound library preparation for virtual screening**

MolScreener takes raw vendor compound libraries — however large, however packaged — and turns them into clean, named, filtered, docking-ready inputs for two workflows:

| Output | Use with |
|--------|----------|
| 3D SDF files (one molecule per record, vendor ID on title line) | GNINA · Vina · AutoDock |
| YAML files (receptor + ligand, one file per compound) | Boltz-2 co-folding & affinity prediction |

Every compound keeps its original catalogue ID throughout — in the output filename, the mol-block title line, and a master index TSV — so hits are always traceable back to the vendor.

---

## Features

- **Any vendor format** — ZIP archives, plain SDF, gzip-compressed SDF (`.sdf.gz`), SMILES files (`.smi`), delimited spreadsheets (`.csv` / `.tsv`), or whole directories scanned recursively
- **Robust ID extraction** — reads `Catalog_ID`, `idnumber`, `PUBCHEM_EXT_DATASOURCE_REGID`, `MOLPORTID`, and 10+ other vendor tags; falls back to a stable `{source}_{index}` ID so nothing is ever unnamed
- **14 drug-likeness filters** — MW, LogP, HBD, HBA, TPSA, rotatable bonds, fCsp³, QED, SA Score, heavy atoms, formal charge, rings, aromatic rings, stereocentres, plus PAINS and Brenk structural-alert removal
- **Four ready-to-use presets** — Drug Design (Ro5), Lead-Like (Ro4), Fragment-Like (Ro3), Kinase-Focused
- **Server parallelism** — inputs are chunked and distributed across all chosen CPUs; a single 6 GB SDF uses 32 cores just as efficiently as 32 small files
- **Two modes** — interactive TUI (arrow-key menus, good for first use) or fully non-interactive CLI (good for SLURM / HPC job scripts)

---

## Installation

### 1 · Python dependencies

```bash
pip install rdkit pyyaml rich questionary
```

Or with conda:

```bash
conda install -c conda-forge rdkit pyyaml rich questionary
```

### 2 · OpenBabel  *(GNINA mode only — for pH-aware 3D conformer generation)*

```bash
# Linux
sudo apt install openbabel          # Debian / Ubuntu
conda install -c conda-forge openbabel

# macOS
brew install open-babel
```

If OpenBabel is not installed, MolScreener falls back to RDKit ETKDG for 3D generation (no pH protonation).

### 3 · SA Score  *(optional — enables synthetic-accessibility filter)*

The SA Score module ships with RDKit but needs to be on the path.  If you installed RDKit via conda or pip it is usually already available automatically.

---

## Quick start

### Interactive TUI

```bash
python MolScreener.py
```

Five guided steps walk you through everything. Recommended for first use.

### Non-interactive CLI

```bash
# GNINA docking prep — mixed input formats, 32 CPUs
python MolScreener.py \
  --lib /data/Enamine.sdf.gz  Enamine \
  --lib /data/MolPort.zip     MolPort \
  --lib /data/custom.smi      Custom  \
  --mode gnina --cpus 32 --ph 7.4 \
  --output /scratch/gnina_ready

# Boltz-2 affinity prep — whole directory, 16 CPUs
python MolScreener.py \
  --lib /data/library/        MyLib \
  --mode boltz2 --cpus 16 \
  --sequence MTEYKLVVVGAGG... \
  --msa /data/GHRHR_UNPAIRED_MSA.a3m \
  --output /scratch/boltz2_ready
```

### Re-run a saved configuration

The TUI offers to save your settings to JSON at the end of setup. That file can then be used for exact re-runs:

```bash
python MolScreener.py --config molscreener_config.json
```

This is useful for SLURM scripts, repeating a campaign with a different library, or sharing a run configuration with a collaborator.

---

## Accepted input formats

| Format | Extensions | Notes |
|--------|------------|-------|
| ZIP archive | `.zip` | Can contain SDF, SMILES, or CSV files; nested folders are fine |
| Plain SDF | `.sdf` | Standard 2D SDF from any vendor |
| Compressed SDF | `.sdf.gz`, `.gz` | How Enamine and MolPort ship large catalogs |
| SMILES file | `.smi`, `.smiles`, `.ism` | One compound per line; with or without ID column |
| Delimited spreadsheet | `.csv`, `.tsv`, `.txt` | Auto-detects SMILES and ID columns from header |
| Directory | — | Scanned recursively for any of the above |

---

## Vendor ID handling

Catalogue IDs are extracted from SDF properties in priority order:

```
Catalog_ID · idnumber · PUBCHEM_EXT_DATASOURCE_REGID · MOLPORTID · MolPort_ID ·
molport_id · CatalogID · Catalogue_ID · catalogue_id · Vendor_ID · Mcule_ID ·
Compound_ID · CompoundID · REGID · Name · ID · CompoundName · mol-block title line
```

If none of these are present, a stable `{source}_{index}` ID is assigned (e.g. `enamine_sdf_4521`). The `id_source` column in the catalogue index records whether each ID came from a vendor tag (`vendor`) or was generated (`fallback`), so coverage can be audited at a glance.

---

## Filters

| Filter | What it measures |
|--------|-----------------|
| MW | Molecular weight (Da) |
| LogP | Lipophilicity |
| HBD / HBA | H-bond donors / acceptors |
| TPSA | Topological polar surface area — governs membrane permeability |
| Rotatable bonds | Molecular flexibility |
| fCsp³ | Fraction of sp³ carbons — correlates with solubility |
| QED | Overall drug-likeness score (0–1) |
| SA Score | Synthetic accessibility (1–10; lower = easier to make) |
| Heavy atoms | Proxy for molecular size |
| Formal charge | Net charge of the molecule |
| Ring count | Total rings |
| Aromatic rings | Flat π-stacking rings |
| Stereocentres | Chiral centres |
| PAINS | Pan-assay interference — compounds that give false positives across many assays |
| Brenk alerts | Reactive or generally undesirable chemical groups |

### Presets

| Preset | MW | LogP | HBD | HBA | Typical use |
|--------|----|------|-----|-----|-------------|
| Drug Design (Ro5) | ≤ 500 | ≤ 5 | ≤ 5 | ≤ 10 | Oral drug candidates |
| Lead-Like (Ro4) | ≤ 400 | ≤ 4 | ≤ 4 | ≤ 8 | Optimisation hits |
| Fragment-Like (Ro3) | ≤ 300 | ≤ 3 | ≤ 3 | ≤ 6 | Fragment-based drug discovery |
| Kinase-Focused | 300–550 | ≤ 5.5 | ≤ 5 | ≤ 10 | ATP-binding pocket compounds |

---

## Outputs

### GNINA mode

```
output_dir/
├── Enamine_gnina_part1.sdf      ← 3D SDF, vendor ID on title line + <CATALOGUE_ID> tag
├── Enamine_gnina_part2.sdf
├── Enamine_catalogue_index.tsv  ← catalogue_id · id_source · smiles · sdf_file
├── MolPort_gnina_part1.sdf
└── MolPort_catalogue_index.tsv
```

### Boltz-2 mode

```
output_dir/
├── Enamine_boltz2_inputs/
│   ├── Z1234567890.yaml         ← one YAML per compound, named by catalogue ID
│   ├── Z1234567891.yaml
│   └── ...
├── Enamine_catalogue_index.tsv  ← catalogue_id · id_source · smiles · yaml_file
└── ...
```

Each Boltz-2 YAML follows the canonical schema:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MTEYKLVVV...
      msa: /path/to/RECEPTOR_UNPAIRED_MSA.a3m
  - ligand:
      id: B
      smiles: 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'
properties:
  - affinity:
      binder: B
      target: A
```

---

## CLI reference

```
python MolScreener.py [options]

Input
  --lib PATH LABEL       Library path (ZIP/SDF/SMILES/CSV/directory) + short label.
                         Repeatable — add as many libraries as needed.
  --config JSON          Path to a saved JSON config. Skips all other flags.

Mode
  --mode gnina|boltz2    Output format (default: gnina)

GNINA options
  --ph FLOAT             Protonation pH for OpenBabel (default: 7.0)
  --batch INT            Max molecules per output SDF file (default: 100000)

Boltz-2 options
  --sequence SEQ         Receptor amino-acid sequence (required for boltz2)
  --msa FILE             Path to receptor MSA .a3m file — must be a specific file,
                         not a directory; unpaired MSA is correct for a single chain
  --protein-id CHAR      Receptor chain ID in YAML (default: A)
  --ligand-id CHAR       Ligand chain ID in YAML (default: B)

Compute
  --cpus INT             Worker processes to use (default: all CPUs minus one)
  --chunk INT            Molecules per parallel task (default: 2000)

Output
  --output DIR           Output directory (created if it does not exist)
```

---

## Example configuration file

See [`examples/example_config.json`](examples/example_config.json) for a complete annotated configuration and [`examples/test_library.smi`](examples/test_library.smi) for a small test library you can run immediately.

---

## Requirements

- Python ≥ 3.10
- RDKit
- PyYAML
- rich
- questionary
- OpenBabel *(optional, recommended for GNINA mode)*

---

## License

MIT — see [LICENSE](LICENSE).

---

## Citation

If you use MolScreener in published work, please cite the tool and the underlying software it depends on:

- **RDKit**: https://www.rdkit.org
- **OpenBabel**: O'Boyle *et al.*, J. Cheminform. 2011, 3, 33
- **Boltz-2**: Wohlwend *et al.*, 2024 (if using Boltz-2 mode)
- **GNINA**: McNutt *et al.*, J. Cheminform. 2021, 13, 43 (if using GNINA mode)
