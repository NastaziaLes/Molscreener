#!/usr/bin/env python3
"""
MolScreener — Compound Library Preparation for Virtual Screening
════════════════════════════════════════════════════════════════════

Takes a raw vendor compound library (however large, however packaged) and
turns it into clean, named, filtered, docking-ready inputs in two modes:

  GNINA / Vina / AutoDock   →  3D SDF files, one molecule per record,
                               vendor ID stamped on the title line and tag

  Boltz-2 co-folding        →  one YAML per compound (receptor + ligand),
                               ready to hand to  boltz predict  as a folder

Run it interactively (recommended for first use):
  python MolScreener.py

Or non-interactively with a saved config (good for job-schedulers, SLURM, etc.):
  python MolScreener.py --config run.json

Or fully inline:
  python MolScreener.py --lib Enamine.sdf.gz Enamine --mode gnina --cpus 32 ...

─────────────────────────────────────────────────────────────────────
WHAT IT ACCEPTS

  .zip          ZIP archive containing SDF / SMILES / CSV files (nested OK)
  .sdf          a plain 2D SDF from a vendor
  .sdf.gz       gzip-compressed SDF  (how large libraries are normally shipped)
  .smi / .smiles / .ism   SMILES file, one compound per line
  .csv / .tsv / .txt      delimited file with a SMILES column + an ID column
  a directory   any mix of the above, scanned recursively

VENDOR ID HANDLING

  Every output — SDF title line, YAML filename, TSV index — carries the
  compound's original catalogue ID.  IDs are read from vendor SDF tags in
  priority order:

    Catalog_ID · idnumber · PUBCHEM_EXT_DATASOURCE_REGID · MOLPORTID ·
    molport_id · CatalogID · Vendor_ID · Mcule_ID · ID · CompoundName ·
    mol-block title line (last resort)

  If none of those are present, a stable  {source}_{index}  ID is assigned
  automatically — nothing is ever unnamed, and the index TSV records whether
  each ID came from a vendor tag ("vendor") or was generated ("fallback").

PARALLELISM

  Inputs are split into fixed-size molecule chunks and spread across the CPUs
  you choose.  A single 6 GB SDF uses all 32 cores just as well as 32 small
  files.  Memory stays flat regardless of library size.

OUTPUTS

  GNINA mode
    {label}_gnina_part1.sdf, part2.sdf, …   (batched 3D SDF, one ID per title)
    {label}_catalogue_index.tsv              (catalogue_id · id_source · smiles · file)

  Boltz-2 mode
    {label}_boltz2_inputs/
      {catalogue_id}.yaml                    (one file per compound, named by ID)
    {label}_catalogue_index.tsv              (catalogue_id · id_source · smiles · yaml)
"""

from __future__ import annotations
import argparse, concurrent.futures, csv, gzip, io, json, multiprocessing
import os, re, shutil, subprocess, sys, tempfile, time, zipfile
from dataclasses import dataclass, asdict, field
from enum import Enum
from pathlib import Path
from typing import Callable, Iterable, Iterator, TextIO

# ── optional deps ────────────────────────────────────────────────────
try:
    import yaml
except ImportError:
    sys.exit("pyyaml missing — pip install pyyaml")

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich import box
    _RICH = True
except ImportError:
    _RICH = False

try:
    import questionary
    from questionary import Style as QStyle
    _Q = True
except ImportError:
    _Q = False

try:
    from rdkit import Chem, RDConfig, RDLogger
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors, rdmolops
    from rdkit.Chem import FilterCatalog, QED as rdQED
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    sys.exit("RDKit missing — conda install -c conda-forge rdkit  (or pip install rdkit)")

try:
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer as _sascorer
    _SA_OK = True
except Exception:
    _SA_OK = False

OPENBABEL  = shutil.which("obabel") or "obabel"
TOTAL_CPUS = multiprocessing.cpu_count()
console    = Console() if _RICH else None

# Recognised input extensions
_SDF_EXT = {".sdf"}
_SMI_EXT = {".smi", ".smiles", ".ism"}
_TAB_EXT = {".csv", ".tsv", ".txt"}
_ALL_MOL_EXT = _SDF_EXT | _SMI_EXT | _TAB_EXT


# ══════════════════════════════════════════════════════════════════════
# §1  DATA STRUCTURES
# ══════════════════════════════════════════════════════════════════════

class OutputMode(str, Enum):
    gnina  = "gnina"
    boltz2 = "boltz2"

class PresetName(str, Enum):
    drug_design    = "drug_design"
    lead_like      = "lead_like"
    fragment_like  = "fragment_like"
    kinase_focused = "kinase_focused"

@dataclass
class ScreeningCriteria:
    min_mw: float = 150;      max_mw: float = 500
    min_logp: float = -2;     max_logp: float = 5
    min_hbd: int = 0;         max_hbd: int = 5
    min_hba: int = 0;         max_hba: int = 10
    min_tpsa: float = 0;      max_tpsa: float = 140
    min_rot_bonds: int = 0;   max_rot_bonds: int = 10
    min_fcsp3: float = 0.2;   max_fcsp3: float = 1.0
    min_qed: float = 0.4;     max_qed: float = 1.0
    min_sa: float = 1.0;      max_sa: float = 6.0
    min_heavy: int = 10;      max_heavy: int = 50
    min_charge: int = -2;     max_charge: int = 2
    min_rings: int = 1;       max_rings: int = 5
    min_arom: int = 0;        max_arom: int = 3
    min_stereo: int = 0;      max_stereo: int = 4
    remove_pains: bool = True
    remove_brenk: bool = True

@dataclass
class Boltz2Config:
    protein_sequence: str = ""
    msa_file:         str = ""
    protein_id:       str = "A"
    ligand_id:        str = "B"

@dataclass
class LibraryEntry:
    label: str
    path:  str   # ZIP, SDF, SMILES, CSV/TSV, .gz, or a directory

@dataclass
class RunConfig:
    libraries:  list[LibraryEntry] = field(default_factory=list)
    criteria:   ScreeningCriteria  = field(default_factory=ScreeningCriteria)
    mode:       OutputMode         = OutputMode.gnina
    output_dir: str                = ""
    n_cpus:     int                = max(1, TOTAL_CPUS - 1)
    chunk_size: int                = 500     # molecules per parallel task
    # GNINA options
    ph:         float              = 7.0
    batch_size: int                = 100_000  # molecules per output SDF part
    # Boltz-2 options
    boltz2:     Boltz2Config       = field(default_factory=Boltz2Config)


# ══════════════════════════════════════════════════════════════════════
# §2  PRESETS
# ══════════════════════════════════════════════════════════════════════

def apply_preset(name: PresetName) -> ScreeningCriteria:
    if name == PresetName.drug_design:
        return ScreeningCriteria()  # dataclass defaults == Ro5 drug design
    elif name == PresetName.lead_like:
        return ScreeningCriteria(
            min_mw=150,    max_mw=400,    min_logp=-2,  max_logp=4,
            min_hbd=0,     max_hbd=4,     min_hba=0,    max_hba=8,
            min_tpsa=0,    max_tpsa=120,  min_rot_bonds=0, max_rot_bonds=7,
            min_fcsp3=0.2, max_fcsp3=1.0, min_qed=0.5,  max_qed=1.0,
            min_sa=1.0,    max_sa=5.0,    min_heavy=10, max_heavy=40,
            min_charge=-1, max_charge=1,  min_rings=1,  max_rings=4,
            min_arom=0,    max_arom=3,    min_stereo=0, max_stereo=3,
            remove_pains=True, remove_brenk=True)
    elif name == PresetName.fragment_like:
        return ScreeningCriteria(
            min_mw=100,    max_mw=300,    min_logp=-2,  max_logp=3,
            min_hbd=0,     max_hbd=3,     min_hba=0,    max_hba=6,
            min_tpsa=0,    max_tpsa=90,   min_rot_bonds=0, max_rot_bonds=5,
            min_fcsp3=0.1, max_fcsp3=1.0, min_qed=0.3,  max_qed=1.0,
            min_sa=1.0,    max_sa=4.0,    min_heavy=7,  max_heavy=27,
            min_charge=-1, max_charge=1,  min_rings=1,  max_rings=3,
            min_arom=0,    max_arom=2,    min_stereo=0, max_stereo=2,
            remove_pains=True, remove_brenk=True)
    else:  # kinase_focused
        return ScreeningCriteria(
            min_mw=300,    max_mw=550,    min_logp=1.0, max_logp=5.5,
            min_hbd=0,     max_hbd=5,     min_hba=2,    max_hba=10,
            min_tpsa=40,   max_tpsa=130,  min_rot_bonds=2, max_rot_bonds=10,
            min_fcsp3=0.1, max_fcsp3=0.7, min_qed=0.35, max_qed=1.0,
            min_sa=1.0,    max_sa=6.0,    min_heavy=22, max_heavy=45,
            min_charge=-1, max_charge=1,  min_rings=2,  max_rings=6,
            min_arom=1,    max_arom=4,    min_stereo=0, max_stereo=3,
            remove_pains=True, remove_brenk=True)


# ══════════════════════════════════════════════════════════════════════
# §3  VENDOR-ID EXTRACTION  (the heart of "names always printed")
# ══════════════════════════════════════════════════════════════════════

# Priority order. Enamine ships `idnumber` / `Catalog_ID`; MolPort ships
# `MOLPORTID` / `PUBCHEM_EXT_DATASOURCE_REGID`. We try the most specific
# vendor tags first and only fall back to generic ones.
_ID_TAGS = [
    "Catalog_ID",                     # Enamine
    "idnumber",                       # Enamine (alt)
    "PUBCHEM_EXT_DATASOURCE_REGID",   # MolPort
    "MOLPORTID",                      # MolPort
    "MolPort_ID",
    "molport_id",
    "CatalogID",
    "Catalogue_ID",
    "catalogue_id",
    "Vendor_ID",
    "Mcule_ID",
    "Compound_ID",
    "CompoundID",
    "REGID",
    "Name",
    "ID",
    "CompoundName",
]

# Title lines that are program stamps / junk, never real IDs.
_TITLE_JUNK = re.compile(r'^(Mrv\d|.*RDKit.*|.*OpenBabel.*|.*Marvin.*|\d+$|\s*$)',
                         re.IGNORECASE)


def extract_vendor_id(mol) -> tuple[str, bool]:
    """
    Return (id, is_vendor).

    is_vendor=True  → a real catalogue ID was found in a tag or the title.
    is_vendor=False → caller must apply a deterministic positional fallback;
                      this function returns ("", False) in that case.

    Never raises, never returns a non-empty junk title.
    """
    for tag in _ID_TAGS:
        if mol.HasProp(tag):
            val = mol.GetProp(tag).strip()
            if val:
                return val, True
    if mol.HasProp("_Name"):
        title = mol.GetProp("_Name").strip()
        if title and not _TITLE_JUNK.match(title):
            return title, True
    return "", False


def _make_catalogs(remove_pains: bool, remove_brenk: bool):
    p_cat = b_cat = None
    if remove_pains:
        p = FilterCatalogParams()
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        p_cat = FilterCatalog.FilterCatalog(p)
    if remove_brenk:
        p = FilterCatalogParams()
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        b_cat = FilterCatalog.FilterCatalog(p)
    return p_cat, b_cat


def _passes(mol, c: ScreeningCriteria, pains_cat, brenk_cat) -> bool:
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rot  = Descriptors.NumRotatableBonds(mol)
    fsp3 = Descriptors.FractionCSP3(mol)
    if not (c.min_mw   <= mw   <= c.max_mw):    return False
    if not (c.min_logp <= logp <= c.max_logp):  return False
    if not (c.min_hbd  <= hbd  <= c.max_hbd):    return False
    if not (c.min_hba  <= hba  <= c.max_hba):    return False
    if not (c.min_tpsa <= tpsa <= c.max_tpsa):   return False
    if not (c.min_rot_bonds <= rot  <= c.max_rot_bonds): return False
    if not (c.min_fcsp3     <= fsp3 <= c.max_fcsp3):     return False

    heavy  = mol.GetNumHeavyAtoms()
    charge = rdmolops.GetFormalCharge(mol)
    rings  = rdMolDescriptors.CalcNumRings(mol)
    arom   = rdMolDescriptors.CalcNumAromaticRings(mol)
    stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if not (c.min_heavy  <= heavy  <= c.max_heavy):   return False
    if not (c.min_charge <= charge <= c.max_charge):  return False
    if not (c.min_rings  <= rings  <= c.max_rings):   return False
    if not (c.min_arom   <= arom   <= c.max_arom):    return False
    if not (c.min_stereo <= stereo <= c.max_stereo):  return False

    if not (c.min_qed <= rdQED.qed(mol) <= c.max_qed): return False

    if _SA_OK and (c.min_sa > 1.0 or c.max_sa < 10.0):
        try:
            sa = _sascorer.calculateScore(mol)
            if not (c.min_sa <= sa <= c.max_sa): return False
        except Exception:
            pass

    if pains_cat and pains_cat.HasMatch(mol): return False
    if brenk_cat and brenk_cat.HasMatch(mol): return False
    return True


# ══════════════════════════════════════════════════════════════════════
# §4  WORKER  (runs in child processes — must be top-level & picklable)
# ══════════════════════════════════════════════════════════════════════
#
# A "chunk" is a tuple: (fmt, payload, index_base, src_tag)
#   fmt == "sdf" : payload is SDF text holding up to chunk_size molecules
#   fmt == "smi" : payload is a list of (raw_id_or_None, smiles)
#   index_base   : running molecule index for deterministic fallback IDs
#   src_tag      : sanitised source name, used to build fallback IDs
#
# Returns: list of (catalogue_id, id_source, smiles)
#   id_source in {"vendor", "fallback"}

def worker_process_chunk(chunk, criteria: ScreeningCriteria):
    fmt, payload, index_base, src_tag = chunk
    RDLogger.DisableLog('rdApp.*')
    pains_cat, brenk_cat = _make_catalogs(criteria.remove_pains,
                                          criteria.remove_brenk)
    out: list[tuple[str, str, str]] = []

    if fmt == "sdf":
        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".sdf")
        try:
            with os.fdopen(tmp_fd, "w", encoding="utf-8", errors="replace") as fh:
                fh.write(payload)
            suppl = Chem.SDMolSupplier(tmp_path, removeHs=True, sanitize=True)
            for i, mol in enumerate(suppl):
                pos = index_base + i
                if mol is None:
                    continue
                cid, is_vendor = extract_vendor_id(mol)
                if not cid:
                    cid, is_vendor = f"{src_tag}_{pos}", False
                rec = _finalize(mol_from="mol", mol=mol, cid=cid,
                                is_vendor=is_vendor, criteria=criteria,
                                pains_cat=pains_cat, brenk_cat=brenk_cat)
                if rec:
                    out.append(rec)
        finally:
            try: os.unlink(tmp_path)
            except OSError: pass

    else:  # "smi"
        for i, (raw_id, smiles) in enumerate(payload):
            pos = index_base + i
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if mol is None:
                continue
            if raw_id and str(raw_id).strip():
                cid, is_vendor = str(raw_id).strip(), True
            else:
                cid, is_vendor = f"{src_tag}_{pos}", False
            rec = _finalize(mol_from="smiles", mol=mol, cid=cid,
                            is_vendor=is_vendor, criteria=criteria,
                            pains_cat=pains_cat, brenk_cat=brenk_cat)
            if rec:
                out.append(rec)

    return out


def _finalize(*, mol_from, mol, cid, is_vendor, criteria, pains_cat, brenk_cat):
    """Canonicalise → re-parse → filter. Returns (cid, id_source, smiles) or None."""
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None
    if not smiles:
        return None
    if mol_from == "mol":
        mol_c = Chem.MolFromSmiles(smiles)   # consistent descriptors
        if mol_c is None:
            return None
    else:
        mol_c = mol
    if _passes(mol_c, criteria, pains_cat, brenk_cat):
        return (cid, "vendor" if is_vendor else "fallback", smiles)
    return None


# ══════════════════════════════════════════════════════════════════════
# §5  INPUT RESOLUTION  (zip / sdf / gz / smi / csv / directory → shards)
# ══════════════════════════════════════════════════════════════════════

@dataclass
class Shard:
    """One readable molecule source plus how to open it as text."""
    name:   str                       # display name
    fmt:    str                       # "sdf" or "smi"
    opener: Callable[[], TextIO]      # returns a fresh text stream


def _strip_gz(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def _fmt_for(name: str) -> str | None:
    ext = Path(_strip_gz(name)).suffix.lower()
    if ext in _SDF_EXT: return "sdf"
    if ext in _SMI_EXT: return "smi"
    if ext in _TAB_EXT: return "smi"   # delimited → parsed as smi rows
    return None


def _open_plain(path: str) -> TextIO:
    if path.lower().endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _open_zip_member(zip_path: str, member: str) -> TextIO:
    arc = zipfile.ZipFile(zip_path, "r")
    raw = arc.open(member, "r")
    if member.lower().endswith(".gz"):
        raw = gzip.open(raw, "rb")
    wrapper = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    # Keep the archive alive for the lifetime of the wrapper.
    wrapper._arc = arc  # type: ignore[attr-defined]
    return wrapper


def resolve_shards(path: str) -> list[Shard]:
    """Expand a library path into a flat list of readable Shards."""
    p = Path(path)
    shards: list[Shard] = []

    if p.is_dir():
        for f in sorted(p.rglob("*")):
            if f.is_file():
                shards.extend(resolve_shards(str(f)))
        return shards

    low = p.name.lower()

    if low.endswith(".zip"):
        with zipfile.ZipFile(path, "r") as arc:
            for info in arc.infolist():
                if info.is_dir() or info.file_size == 0:
                    continue
                fmt = _fmt_for(info.filename)
                if fmt is None:
                    continue
                member = info.filename
                shards.append(Shard(
                    name=f"{p.name}:{member}",
                    fmt=fmt,
                    opener=(lambda zp=path, m=member: _open_zip_member(zp, m)),
                ))
        return shards

    fmt = _fmt_for(low)
    if fmt is not None:
        shards.append(Shard(name=p.name, fmt=fmt,
                            opener=(lambda fp=path: _open_plain(fp))))
    return shards


def src_tag_of(name: str) -> str:
    """Filesystem-safe tag used to build deterministic fallback IDs."""
    base = Path(name.replace(":", "_")).name
    base = re.sub(r'\.(sdf|smi|smiles|ism|csv|tsv|txt|gz)$', '', base,
                  flags=re.IGNORECASE)
    return re.sub(r'[^A-Za-z0-9\-_]', '_', base) or "lib"


# ── streaming splitters ──────────────────────────────────────────────

def _iter_sdf_chunks(stream: TextIO, chunk_size: int):
    """Yield (sdf_text, index_base) without loading the whole file."""
    buf: list[str] = []
    count = 0
    base = 0
    for line in stream:
        buf.append(line)
        if line.startswith("$$$$"):
            count += 1
            if count >= chunk_size:
                yield "".join(buf), base
                base += count
                buf, count = [], 0
    if any(s.strip() for s in buf):
        yield "".join(buf), base


def _sniff_smiles_rows(stream: TextIO) -> Iterator[tuple[str | None, str]]:
    """
    Parse a .smi / .csv / .tsv / .txt file into (id, smiles) rows.

    Handles: headerless "SMILES <id>" (space/tab), and delimited files with a
    header naming a SMILES column and an ID column.
    """
    first = stream.readline()
    if not first:
        return
    delim = "\t" if "\t" in first else ("," if "," in first else None)

    def split(line: str) -> list[str]:
        if delim:
            return next(csv.reader([line], delimiter=delim))
        return line.split()

    cols0 = split(first.rstrip("\n"))
    # Is the first row a header (no valid SMILES in any cell)?
    header = None
    smi_idx, id_idx = 0, 1
    if cols0 and Chem.MolFromSmiles(cols0[0]) is None:
        looks_like_header = any(re.search(r'smiles|smile|structure', c, re.I)
                                for c in cols0)
        if looks_like_header or all(Chem.MolFromSmiles(c) is None for c in cols0):
            header = [c.strip().lower() for c in cols0]
            smi_idx = next((i for i, c in enumerate(header)
                            if re.search(r'smiles|smile|structure', c)), 0)
            id_idx = next((i for i, c in enumerate(header)
                           if re.search(r'(catalog|id|name|compound)', c)
                           and i != smi_idx), None)
    if header is None:
        # first row is data
        smi = cols0[smi_idx] if len(cols0) > smi_idx else ""
        cid = cols0[id_idx] if (id_idx is not None and len(cols0) > id_idx) else None
        if smi:
            yield (cid, smi)

    for line in stream:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        cols = split(line)
        if not cols or len(cols) <= smi_idx:
            continue
        smi = cols[smi_idx].strip()
        cid = (cols[id_idx].strip()
               if (id_idx is not None and len(cols) > id_idx) else None)
        if smi:
            yield (cid, smi)


def _iter_smi_chunks(rows: Iterable[tuple[str | None, str]], chunk_size: int):
    chunk: list[tuple[str | None, str]] = []
    base = 0
    for row in rows:
        chunk.append(row)
        if len(chunk) >= chunk_size:
            yield chunk, base
            base += len(chunk)
            chunk = []
    if chunk:
        yield chunk, base


def enumerate_chunks(library: LibraryEntry, chunk_size: int):
    """Yield ready-to-dispatch chunks for every shard in a library."""
    for shard in resolve_shards(library.path):
        tag = src_tag_of(shard.name)
        stream = shard.opener()
        try:
            if shard.fmt == "sdf":
                for text, base in _iter_sdf_chunks(stream, chunk_size):
                    yield ("sdf", text, base, tag)
            else:
                rows = _sniff_smiles_rows(stream)
                for rows_chunk, base in _iter_smi_chunks(rows, chunk_size):
                    yield ("smi", rows_chunk, base, tag)
        finally:
            try: stream.close()
            except Exception: pass


# ══════════════════════════════════════════════════════════════════════
# §6  PARALLEL DISPATCH  (bounded memory, all CPUs busy)
# ══════════════════════════════════════════════════════════════════════

def dispatch(chunks: Iterator, n_cpus: int, criteria: ScreeningCriteria,
             stop_flag: list[bool], progress_cb=None):
    """
    Stream chunks to a process pool with a read-ahead buffer large enough
    to keep all CPUs busy even when input comes from a single large file.

    The read-ahead window is n_cpus * 4 — so with 20 CPUs, up to 80 chunks
    are submitted at once.  Workers pick them up as they finish; the main
    process keeps reading ahead to refill the queue.  Memory stays bounded
    because each chunk holds at most chunk_size molecules.
    """
    # 4× gives enough runway for the file reader to stay ahead of the workers.
    # With a single large SDF and 20 CPUs, 2× was too small — workers finished
    # chunks faster than the reader could produce the next one.
    max_inflight = max(4, n_cpus * 4)
    it = iter(chunks)
    done = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cpus) as pool:
        inflight = set()
        # Pre-fill the entire window before waiting for any result
        for _ in range(max_inflight):
            try:
                inflight.add(pool.submit(worker_process_chunk, next(it), criteria))
            except StopIteration:
                break

        while inflight:
            completed, inflight = concurrent.futures.wait(
                inflight, return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in completed:
                try:
                    for rec in fut.result():
                        yield rec
                except Exception as ex:
                    print(f"  [!] Worker error: {ex}", file=sys.stderr)
                done += 1
                if progress_cb:
                    progress_cb(done)
            if stop_flag[0]:
                for f in inflight:
                    f.cancel()
                break
            # Refill: for every slot that freed up, read ahead one more chunk
            for _ in range(len(completed)):
                try:
                    inflight.add(pool.submit(worker_process_chunk, next(it), criteria))
                except StopIteration:
                    break


# ══════════════════════════════════════════════════════════════════════
# §7  OUTPUT WRITERS
# ══════════════════════════════════════════════════════════════════════

def _stamp_sdf_id(sdf_block: str, catalogue_id: str) -> str:
    lines = sdf_block.splitlines(keepends=True)
    if lines:
        lines[0] = catalogue_id + "\n"
    text = "".join(lines)
    text = re.sub(r'>\s*<CATALOGUE_ID>\s*\n.*?\n\n', '', text, flags=re.DOTALL)
    tag = f"> <CATALOGUE_ID>\n{catalogue_id}\n\n"
    if "$$$$" in text:
        text = text.replace("$$$$", tag + "$$$$", 1)
    else:
        text = text.rstrip() + "\n" + tag + "$$$$\n"
    return text


def smiles_to_3d_sdf(smiles: str, catalogue_id: str,
                     ph: float, use_ob: bool) -> str | None:
    if use_ob:
        try:
            proc = subprocess.run(
                [OPENBABEL, "-ismi", "-", "-osdf",
                 "--gen3d", "--ph", str(ph), "--title", catalogue_id],
                input=smiles.encode(), capture_output=True, timeout=90)
            if proc.returncode == 0 and proc.stdout:
                raw = proc.stdout.decode("utf-8", errors="replace")
                return _stamp_sdf_id(raw, catalogue_id)
        except Exception:
            pass
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol.SetProp("_Name", catalogue_id)
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(mol)
        block = Chem.MolToMolBlock(mol)
        return block + f"> <CATALOGUE_ID>\n{catalogue_id}\n\n$$$$\n"
    except Exception:
        return None


def make_boltz2_yaml(smiles: str, cfg: Boltz2Config) -> str:
    """
    Build a schema-clean Boltz-2 YAML for one ligand.

    The vendor catalogue ID is intentionally NOT written into the body — it is
    not part of the Boltz schema and Boltz never reads it.  Traceability is
    kept via the YAML filename (= sanitised catalogue ID) and the master
    *_catalogue_index.tsv, which holds the original, unmodified ID.
    """
    doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": cfg.protein_id,
                         "sequence": cfg.protein_sequence,
                         "msa": cfg.msa_file or f"{cfg.protein_id}_MSA.a3m"}},
            {"ligand": {"id": cfg.ligand_id,
                        "smiles": smiles}},
        ],
        "properties": [{"affinity": {"binder": cfg.ligand_id,
                                     "target": cfg.protein_id}}],
    }
    return yaml.dump(doc, allow_unicode=True, sort_keys=False,
                     default_flow_style=False)


def write_gnina(stream, out_dir, lib_label, ph, batch_size,
                stop_flag, progress_cb=None) -> tuple[int, int, int]:
    """Returns (saved, failed_3d, fallback_ids)."""
    use_ob = bool(shutil.which(OPENBABEL))
    if not use_ob:
        _log("warn", "obabel not found — using RDKit ETKDG (no pH correction)")

    tsv_path = os.path.join(out_dir, f"{lib_label}_catalogue_index.tsv")
    tsv = open(tsv_path, "w", encoding="utf-8")
    tsv.write("catalogue_id\tid_source\tsmiles\tsdf_file\n")

    saved = failed = fallback = 0
    part, in_part = 1, 0
    part_name = f"{lib_label}_gnina_part{part}.sdf"
    fh = open(os.path.join(out_dir, part_name), "w", encoding="utf-8")
    try:
        for cid, id_source, smiles in stream:
            if stop_flag[0]:
                break
            if id_source == "fallback":
                fallback += 1
            block = smiles_to_3d_sdf(smiles, cid, ph, use_ob)
            if block is None:
                tsv.write(f"{cid}\t{id_source}\t{smiles}\tFAILED_3D_GEN\n")
                failed += 1
                continue
            fh.write(block)
            tsv.write(f"{cid}\t{id_source}\t{smiles}\t{part_name}\n")
            saved += 1
            in_part += 1
            if in_part >= batch_size:
                fh.flush(); fh.close()
                part += 1; in_part = 0
                part_name = f"{lib_label}_gnina_part{part}.sdf"
                fh = open(os.path.join(out_dir, part_name), "w", encoding="utf-8")
            if progress_cb:
                progress_cb(saved, failed, cid)
    finally:
        fh.flush(); fh.close()
        tsv.flush(); tsv.close()
        last = os.path.join(out_dir, part_name)
        if in_part == 0 and part > 1 and os.path.exists(last):
            os.remove(last)
    return saved, failed, fallback


def write_boltz2(stream, out_dir, lib_label, cfg,
                 stop_flag, progress_cb=None) -> tuple[int, int]:
    """Returns (saved, fallback_ids)."""
    yaml_dir = os.path.join(out_dir, f"{lib_label}_boltz2_inputs")
    os.makedirs(yaml_dir, exist_ok=True)
    tsv_path = os.path.join(out_dir, f"{lib_label}_catalogue_index.tsv")
    tsv = open(tsv_path, "w", encoding="utf-8")
    tsv.write("catalogue_id\tid_source\tsmiles\tyaml_file\n")

    saved = fallback = 0
    seen: dict[str, int] = {}
    try:
        for cid, id_source, smiles in stream:
            if stop_flag[0]:
                break
            if id_source == "fallback":
                fallback += 1
            safe = re.sub(r'[^A-Za-z0-9\-_]', '_', cid)
            # guard against filename collisions across sources
            if safe in seen:
                seen[safe] += 1
                safe = f"{safe}__{seen[safe]}"
            else:
                seen[safe] = 0
            yaml_filename = f"{safe}.yaml"
            with open(os.path.join(yaml_dir, yaml_filename), "w",
                      encoding="utf-8") as f:
                f.write(make_boltz2_yaml(smiles, cfg))
            tsv.write(f"{cid}\t{id_source}\t{smiles}\t{yaml_filename}\n")
            saved += 1
            if saved % 500 == 0:
                tsv.flush()
            if progress_cb:
                progress_cb(saved, cid)
    finally:
        tsv.flush(); tsv.close()
    return saved, fallback


# ══════════════════════════════════════════════════════════════════════
# §8  CORE PIPELINE
# ══════════════════════════════════════════════════════════════════════

def validate_msa_path(msa: str) -> None:
    """Reject a directory or missing file before a (possibly long) Boltz run."""
    if not msa or msa == "empty":
        return
    if os.path.isdir(msa):
        sys.exit(
            f"MSA path is a directory: {msa}\n"
            "Boltz needs a single .a3m FILE for the receptor, not a folder. "
            "Point --msa (or the config 'msa') at the specific file, "
            "e.g. .../GHRHR_UNPAIRED_MSA.a3m")
    if not os.path.isfile(msa):
        sys.exit(f"MSA file not found: {msa}")


def run_pipeline(cfg: RunConfig, stop_flag: list[bool]) -> None:
    if cfg.mode == OutputMode.boltz2:
        validate_msa_path(cfg.boltz2.msa_file)
    os.makedirs(cfg.output_dir, exist_ok=True)

    if console:
        console.print(Panel(
            f"[bold]Mode[/bold]       {cfg.mode.value.upper()}   "
            f"[dim]({'2D → 3D SDF for docking' if cfg.mode == OutputMode.gnina else '2D → YAML for co-folding'})[/dim]\n"
            f"[bold]Libraries[/bold]  {len(cfg.libraries)}\n"
            f"[bold]CPUs[/bold]       {cfg.n_cpus} / {TOTAL_CPUS}   "
            f"[dim]({cfg.chunk_size:,} molecules per parallel task)[/dim]\n"
            f"[bold]Output[/bold]     {cfg.output_dir}",
            title="[bold cyan]▶  Running[/bold cyan]",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))
    else:
        _log("info", f"Running {cfg.mode.value.upper()} · {len(cfg.libraries)} "
                     f"library(ies) · {cfg.n_cpus}/{TOTAL_CPUS} CPUs · "
                     f"output → {cfg.output_dir}")

    stats: list[dict] = []
    t0 = time.time()

    for lib_idx, lib in enumerate(cfg.libraries, 1):
        prefix = f"[{lib_idx}/{len(cfg.libraries)}] {lib.label}"
        rec = {"label": lib.label, "files": 0, "kept": 0,
               "failed": 0, "fallback": 0}
        _log("step", f"{prefix} — reading input…")

        try:
            shards = resolve_shards(lib.path)
        except Exception as ex:
            _log("error", f"{prefix}: cannot read input — {ex}")
            stats.append(rec); continue
        if not shards:
            _log("warn", f"{prefix}: no .sdf / .smi / .csv content found.")
            stats.append(rec); continue

        rec["files"] = len(shards)
        kinds = ", ".join(sorted({s.fmt for s in shards}))
        _log("info", f"{prefix}: {len(shards)} source file(s) [{kinds}] — filtering…")

        def filter_cb(done, lbl=prefix):
            _progress(f"{lbl}: filtered {done} chunk(s)…")

        chunks = enumerate_chunks(lib, cfg.chunk_size)
        stream = dispatch(chunks, cfg.n_cpus, cfg.criteria, stop_flag, filter_cb)

        if cfg.mode == OutputMode.gnina:
            def gnina_cb(saved, failed, cid):
                _progress(f"{prefix}: {saved:,} written"
                          + (f", {failed} failed" if failed else "")
                          + f"   · last ID: {cid}")
            saved, failed, fb = write_gnina(
                stream, cfg.output_dir, lib.label,
                cfg.ph, cfg.batch_size, stop_flag, gnina_cb)
            rec.update(kept=saved, failed=failed, fallback=fb)
            _erase_progress()
            _log("ok", f"{prefix}: {saved:,} molecules → 3D SDF"
                 + (f"  ({failed} failed 3D generation)" if failed else ""))
        else:
            def boltz_cb(saved, cid):
                _progress(f"{prefix}: {saved:,} YAML(s)   · last ID: {cid}")
            saved, fb = write_boltz2(
                stream, cfg.output_dir, lib.label,
                cfg.boltz2, stop_flag, boltz_cb)
            rec.update(kept=saved, fallback=fb)
            _erase_progress()
            _log("ok", f"{prefix}: {saved:,} molecules → Boltz-2 YAML")

        stats.append(rec)

    _print_run_report(cfg, stats, time.time() - t0, stop_flag[0])


def _print_run_report(cfg: RunConfig, stats: list[dict],
                      elapsed: float, aborted: bool) -> None:
    """A clear end-of-run summary: per-library counts + where the output is."""
    total_kept = sum(s["kept"] for s in stats)
    total_failed = sum(s["failed"] for s in stats)
    total_fb = sum(s["fallback"] for s in stats)
    rate = (total_kept / elapsed) if elapsed > 0 else 0
    out_word = "3D SDF" if cfg.mode == OutputMode.gnina else "YAML"

    if console:
        t = Table(box=box.SIMPLE_HEAVY, border_style="dim", pad_edge=False,
                  header_style="bold cyan")
        t.add_column("Library", style="white", no_wrap=True)
        t.add_column("Files", justify="right", style="dim")
        t.add_column(f"Prepared ({out_word})", justify="right", style="green")
        if cfg.mode == OutputMode.gnina:
            t.add_column("Failed 3D", justify="right", style="yellow")
        t.add_column("Fallback IDs", justify="right", style="yellow")
        for s in stats:
            row = [s["label"], f"{s['files']:,}", f"{s['kept']:,}"]
            if cfg.mode == OutputMode.gnina:
                row.append(f"{s['failed']:,}" if s["failed"] else "—")
            row.append(f"{s['fallback']:,}" if s["fallback"] else "—")
            t.add_row(*row)
        if len(stats) > 1:
            total_row = ["[bold]TOTAL[/bold]", "",
                         f"[bold]{total_kept:,}[/bold]"]
            if cfg.mode == OutputMode.gnina:
                total_row.append(f"{total_failed:,}" if total_failed else "—")
            total_row.append(f"{total_fb:,}" if total_fb else "—")
            t.add_row(*total_row)
        console.print(t)

        head = "[yellow]⚠ Aborted[/yellow]" if aborted else "[bold green]✓ Done[/bold green]"
        lines = [
            f"{head}   {total_kept:,} molecules prepared in {elapsed:.0f}s "
            f"([dim]{rate:,.0f}/s[/dim])",
            f"[bold]Output:[/bold] {cfg.output_dir}",
            f"[bold]Index:[/bold]  one [cyan]<label>_catalogue_index.tsv[/cyan] per library "
            f"(catalogue_id · id_source · smiles · file)",
        ]
        if total_fb:
            lines.append(
                f"[yellow]Note:[/yellow] {total_fb:,} molecule(s) had no vendor ID "
                f"and were given a stable [cyan]{{source}}_{{index}}[/cyan] ID "
                f"(id_source=fallback). Check those inputs if you expected IDs everywhere.")
        console.print(Panel("\n".join(lines), box=box.ROUNDED,
                            border_style=("yellow" if aborted else "green"),
                            padding=(0, 2)))
    else:
        verb = "Aborted." if aborted else "Done."
        print(f"\n{verb}  {total_kept:,} molecules prepared in {elapsed:.0f}s "
              f"({rate:,.0f}/s)")
        for s in stats:
            extra = (f", {s['failed']} failed 3D" if s["failed"] else "")
            extra += (f", {s['fallback']} fallback IDs" if s["fallback"] else "")
            print(f"  {s['label']}: {s['kept']:,} {out_word}{extra}")
        print(f"  Output: {cfg.output_dir}")
        print(f"  Index : <label>_catalogue_index.tsv per library")
        if total_fb:
            print(f"  Note  : {total_fb} molecule(s) used fallback IDs "
                  f"(id_source=fallback in the index).")


# ══════════════════════════════════════════════════════════════════════
# §9  TERMINAL LOGGING HELPERS
# ══════════════════════════════════════════════════════════════════════

_progress_active = False

def _log(level: str, msg: str) -> None:
    global _progress_active
    if _progress_active:
        print(); _progress_active = False
    # Rich markup version / plain-text fallback
    styled = {
        "info":  ("[cyan]·[/cyan]",    "·"),
        "step":  ("[bold blue]›[/bold blue]", "›"),
        "warn":  ("[yellow]⚠[/yellow]", "!"),
        "error": ("[bold red]✗[/bold red]", "✗"),
        "ok":    ("[bold green]✓[/bold green]", "✓"),
        "done":  ("[bold bright_green]★[/bold bright_green]", "★"),
    }
    rich_icon, plain_icon = styled.get(level, ("[dim]·[/dim]", "·"))
    if console:
        console.print(f"  {rich_icon}  {msg}")
    else:
        print(f"  {plain_icon}  {msg}", flush=True)

def _progress(msg: str) -> None:
    global _progress_active
    cols = shutil.get_terminal_size((80, 24)).columns
    text = msg[:cols - 4]
    print(f"\r  {text:<{cols - 4}}", end="", flush=True)
    _progress_active = True

def _erase_progress() -> None:
    global _progress_active
    if _progress_active:
        cols = shutil.get_terminal_size((80, 24)).columns
        print(f"\r{' ' * cols}\r", end="", flush=True)
        _progress_active = False


# ══════════════════════════════════════════════════════════════════════
# §10  INTERACTIVE TUI
# ══════════════════════════════════════════════════════════════════════

Q_STYLE = None
if _Q:
    Q_STYLE = QStyle([
        ("qmark", "fg:#4EA8DE bold"), ("question", "bold"),
        ("answer", "fg:#4EC9B0 bold"), ("pointer", "fg:#4EA8DE bold"),
        ("highlighted", "fg:#4EA8DE bold"), ("selected", "fg:#4EC9B0"),
        ("instruction", "fg:#888888"),
    ])

def _q(prompt_fn, *args, **kwargs):
    if _Q:
        kwargs.setdefault("style", Q_STYLE)
    return prompt_fn(*args, **kwargs).ask()


def tui_collect_libraries() -> list[LibraryEntry]:
    if console:
        console.print(Panel(
            "Point me at one or more compound libraries. I will scan them,\n"
            "extract every molecule, and read the catalogue ID from the vendor tags\n"
            "so nothing is ever left unnamed.\n\n"
            "[bold]Accepted formats:[/bold]\n"
            "  [cyan].zip[/cyan]               ZIP archive "
            "(can contain SDF, SMILES, or CSV files — nested folders are fine)\n"
            "  [cyan].sdf[/cyan]               plain 2D SDF from a vendor\n"
            "  [cyan].sdf.gz[/cyan]            gzip-compressed SDF "
            "(how Enamine and MolPort ship large catalogs)\n"
            "  [cyan].smi / .smiles[/cyan]     SMILES file, one compound per line\n"
            "  [cyan].csv / .tsv[/cyan]        delimited spreadsheet with a SMILES column\n"
            "  [cyan]a folder[/cyan]           scanned recursively for any of the above\n\n"
            "[dim]You can add as many libraries as you like. "
            "Each gets a short label that appears in the output filenames.[/dim]",
            title="  Step 1 of 5 · Input library  ",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))

    libs = []
    while True:
        path = _q(questionary.path,
                  "Path to library (ZIP / SDF / SDF.gz / SMILES / CSV / folder):")
        if not path:
            break
        path = str(Path(path).expanduser().resolve())
        if not os.path.exists(path):
            _log("error", f"Path not found: {path}"); continue
        try:
            shards = resolve_shards(path)
        except Exception as ex:
            _log("error", f"Cannot read {path}: {ex}"); continue
        if not shards:
            _log("warn", f"No recognisable molecule files found in {path} — skipping.")
            continue
        counts = {}
        for s in shards:
            counts[s.fmt] = counts.get(s.fmt, 0) + 1
        desc = ", ".join(f"{v} {k.upper()} file{'s' if v > 1 else ''}"
                         for k, v in sorted(counts.items()))
        _log("ok", f"Found {len(shards)} source file(s): {desc}")

        default_label = Path(path).stem
        label = (_q(questionary.text,
                    "Short label for this library (used in output filenames):",
                    default=default_label) or default_label).strip()
        libs.append(LibraryEntry(label=label, path=path))

        if not _q(questionary.confirm, "Add another library?", default=False):
            break

    if not libs:
        _log("error", "No valid libraries added. Exiting.")
        sys.exit(1)
    return libs


def tui_collect_criteria() -> ScreeningCriteria:
    if console:
        console.print(Panel(
            "Choose which molecules to keep. Start with a preset that fits your\n"
            "project, then optionally tweak individual cutoffs.\n\n"
            "[bold]Presets:[/bold]\n"
            "  [cyan]Drug Design (Ro5)[/cyan]    Lipinski's Rule of 5. "
            "The standard for oral drugs.\n"
            "                       MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10\n\n"
            "  [cyan]Lead-Like (Ro4)[/cyan]      Smaller, more efficient leads "
            "for optimisation campaigns.\n"
            "                       MW ≤ 400, LogP ≤ 4 — stricter than Ro5\n\n"
            "  [cyan]Fragment-Like (Ro3)[/cyan]  Tiny fragments for FBDD. "
            "High ligand efficiency expected.\n"
            "                       MW ≤ 300, LogP ≤ 3\n\n"
            "  [cyan]Kinase-Focused[/cyan]       Tuned for ATP-binding pockets — "
            "slightly larger and more aromatic.\n"
            "                       MW 300–550, more aromatic rings, higher fCsp³\n\n"
            "[dim]All presets also remove PAINS (pan-assay interference compounds)\n"
            "and Brenk structural alerts (reactive / undesirable groups).[/dim]",
            title="  Step 2 of 5 · Filtering profile  ",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))

    preset_choice = _q(questionary.select, "Which preset fits your campaign?",
        choices=[
            questionary.Choice(
                "Drug Design  (Ro5)  — oral drugs, MW ≤ 500",
                value=PresetName.drug_design),
            questionary.Choice(
                "Lead-Like    (Ro4)  — optimisation hits, MW ≤ 400",
                value=PresetName.lead_like),
            questionary.Choice(
                "Fragment-Like (Ro3) — FBDD fragments, MW ≤ 300",
                value=PresetName.fragment_like),
            questionary.Choice(
                "Kinase-Focused      — ATP-pocket compounds, MW 300–550",
                value=PresetName.kinase_focused),
            questionary.Choice(
                "Custom              — start from Drug Design and edit every cutoff",
                value=None),
        ])
    c = apply_preset(preset_choice) if preset_choice else apply_preset(PresetName.drug_design)

    if not _q(questionary.confirm,
              "Would you like to customise any individual cutoffs?", default=False):
        return c

    if console:
        console.print(
            "\n  [dim]Press ENTER to keep the preset value. "
            "Type a new number to override.[/dim]\n")

    def ask_range(label, hint, amin, amax, is_int=False):
        cast = int if is_int else float
        if console:
            console.print(f"  [bold]{label}[/bold]  [dim]{hint}[/dim]")
        rmin = _q(questionary.text, f"    min  (preset: {getattr(c, amin)}):",
                  default=str(getattr(c, amin)))
        rmax = _q(questionary.text, f"    max  (preset: {getattr(c, amax)}):",
                  default=str(getattr(c, amax)))
        try: setattr(c, amin, cast(rmin))
        except ValueError: pass
        try: setattr(c, amax, cast(rmax))
        except ValueError: pass

    if console:
        console.print("  [cyan]── Drug-likeness ──────────────────────────────[/cyan]")
    ask_range("Molecular weight (Da)",  "heavier molecules cross membranes less well",
              "min_mw", "max_mw")
    ask_range("LogP",   "lipophilicity; too high → poor solubility & toxicity",
              "min_logp", "max_logp")
    ask_range("H-bond donors (HBD)",   "–NH, –OH groups",
              "min_hbd", "max_hbd", True)
    ask_range("H-bond acceptors (HBA)", "N, O atoms that accept H-bonds",
              "min_hba", "max_hba", True)
    ask_range("TPSA (Å²)",  "topological polar surface area — governs membrane permeability",
              "min_tpsa", "max_tpsa")
    ask_range("QED",   "overall drug-likeness score 0–1; 1 = ideal",
              "min_qed", "max_qed")
    ask_range("SA Score",  "synthetic accessibility 1–10; lower = easier to make",
              "min_sa", "max_sa")

    if console:
        console.print("  [cyan]── Shape & complexity ─────────────────────────[/cyan]")
    ask_range("Rotatable bonds",  "flexible bonds; too many → poor oral absorption",
              "min_rot_bonds", "max_rot_bonds", True)
    ask_range("fCsp³",  "fraction of sp³ carbons; higher → better solubility",
              "min_fcsp3", "max_fcsp3")
    ask_range("Heavy atoms",  "proxy for molecular size",
              "min_heavy", "max_heavy", True)
    ask_range("Formal charge",  "overall charge of the molecule",
              "min_charge", "max_charge", True)
    ask_range("Ring count",  "total rings",
              "min_rings", "max_rings", True)
    ask_range("Aromatic rings",  "flat, pi-stacking rings",
              "min_arom", "max_arom", True)
    ask_range("Stereocentres",  "chiral centres; more = harder to synthesise",
              "min_stereo", "max_stereo", True)

    if console:
        console.print("  [cyan]── Structural alerts ──────────────────────────[/cyan]")
    c.remove_pains = _q(questionary.confirm,
        "Remove PAINS?  (pan-assay interference — compounds that give false positives "
        "in many assays)", default=c.remove_pains)
    c.remove_brenk = _q(questionary.confirm,
        "Remove Brenk alerts?  (reactive or generally undesirable chemical groups)",
        default=c.remove_brenk)
    return c


def tui_collect_mode() -> tuple[OutputMode, float, int, Boltz2Config | None]:
    if console:
        console.print(Panel(
            "Choose what to generate from your filtered compounds:\n\n"
            "  [bold]GNINA / Vina / AutoDock[/bold]\n"
            "  Each 2D molecule is converted to a 3D conformer "
            "(OpenBabel with pH protonation, or RDKit ETKDG as fallback)\n"
            "  and written to batched SDF files. The catalogue ID appears on\n"
            "  the mol-block title line so every docking hit is traceable.\n\n"
            "  [bold]Boltz-2  (co-folding + affinity)[/bold]\n"
            "  One YAML file per compound is written, containing the receptor\n"
            "  sequence, receptor MSA, and the ligand SMILES. Hand the whole\n"
            "  folder to  [cyan]boltz predict[/cyan]  to score all compounds at once.\n"
            "  The YAML filename is the catalogue ID, so Boltz's output folders\n"
            "  are automatically named after your compounds.",
            title="  Step 3 of 5 · Output format  ",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))

    mode = _q(questionary.select, "What do you want to produce?", choices=[
        questionary.Choice(
            "3D SDF files  — ready for GNINA / Vina / AutoDock docking",
            value=OutputMode.gnina),
        questionary.Choice(
            "Boltz-2 YAMLs — ready for  boltz predict  (co-folding + affinity)",
            value=OutputMode.boltz2),
    ])

    ph, batch, bcfg = 7.0, 100_000, None

    if mode == OutputMode.gnina:
        if console:
            console.print("\n  [dim]3D conformer settings[/dim]")
        try:
            ph = float(_q(questionary.text,
                "  Protonation pH for OpenBabel  [7.0 = physiological]:",
                default="7.0"))
        except ValueError:
            ph = 7.0
        try:
            batch = int(_q(questionary.text,
                "  Max molecules per output SDF file  "
                "[100000 keeps files manageable]:",
                default="100000"))
        except ValueError:
            batch = 100_000

    else:  # boltz2
        if console:
            console.print("\n  [dim]Receptor & MSA settings (same for every compound)[/dim]")
        bcfg = Boltz2Config()

        bcfg.protein_sequence = _q(questionary.text,
            "  Receptor sequence  (single-letter amino acids, no FASTA header):")

        # Require a specific .a3m FILE, not a directory
        if console:
            console.print(
                "\n  [dim]The MSA (multiple sequence alignment) tells Boltz-2 about "
                "receptor evolution.\n"
                "  Point to the specific  [cyan].a3m[/cyan]  file — unpaired is correct "
                "for a single receptor chain.\n"
                "  Leave blank to skip (reduces accuracy but avoids needing a pre-built MSA).[/dim]")

        a3m_filter = lambda p: os.path.isdir(p) or p.lower().endswith(".a3m")
        raw = _q(questionary.path,
                 "  Receptor MSA file  (.a3m, unpaired, blank to skip):",
                 default="", file_filter=a3m_filter)
        while raw:
            p = Path(raw).expanduser()
            if p.is_dir():
                _log("error",
                     f"That is a folder, not a file: {p}\n"
                     "        Boltz-2 needs a single .a3m file "
                     "(e.g. GHRHR_UNPAIRED_MSA.a3m), not the whole MSA directory.")
            elif not p.is_file():
                _log("error", f"File not found: {p}")
            elif p.suffix.lower() != ".a3m" and not _q(questionary.confirm,
                    f"  {p.name} doesn't end in .a3m — use it anyway?",
                    default=False):
                pass
            else:
                bcfg.msa_file = str(p.resolve())
                _log("ok", f"Receptor MSA: {bcfg.msa_file}")
                break
            raw = _q(questionary.path,
                     "  MSA .a3m file (blank to skip):",
                     default="", file_filter=a3m_filter)

        bcfg.protein_id = _q(questionary.text,
            "  Receptor chain ID in the YAML  [A]:", default="A")
        bcfg.ligand_id  = _q(questionary.text,
            "  Ligand chain ID in the YAML   [B]:", default="B")

    return mode, ph, batch, bcfg


def tui_collect_compute() -> tuple[int, int]:
    if console:
        console.print(Panel(
            f"This server has [bold yellow]{TOTAL_CPUS}[/bold yellow] logical CPUs available.\n\n"
            "MolScreener splits every library into molecule chunks and processes\n"
            "them in parallel — so all the CPUs you allocate stay busy, whether\n"
            "your library is one huge SDF or thousands of small files.\n\n"
            "[bold]CPU count[/bold]   How many parallel worker processes to run.\n"
            "              Leaving one CPU free keeps the server responsive.\n\n"
            "[bold]Chunk size[/bold]  Molecules per task. 2000 is a safe default.\n"
            "              Smaller → work spreads more evenly across CPUs.\n"
            "              Larger  → slightly less overhead per task.",
            title="  Step 4 of 5 · Compute resources  ",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))

    default_cpus = max(1, TOTAL_CPUS - 1)
    try:
        n = int(_q(questionary.text,
                   f"  How many CPUs to use?  (1 – {TOTAL_CPUS}, "
                   f"recommended: {default_cpus}):",
                   default=str(default_cpus)))
        n = max(1, min(n, TOTAL_CPUS))
    except (ValueError, TypeError):
        n = default_cpus
    try:
        chunk = int(_q(questionary.text,
                       "  Molecules per parallel task (chunk size)  [500]:",
                       default="500"))
        chunk = max(100, chunk)
    except (ValueError, TypeError):
        chunk = 2000

    _log("ok", f"Will use {n} CPU(s) with chunks of {chunk:,} molecules.")
    return n, chunk


def tui_collect_output() -> str:
    if console:
        console.print(Panel(
            "Where should the prepared files land?\n\n"
            "MolScreener will write everything here:\n"
            "  · SDF or YAML output files\n"
            "  · a [cyan]<label>_catalogue_index.tsv[/cyan] per library\n"
            "    mapping every catalogue ID → SMILES → output filename\n\n"
            "[dim]The folder is created automatically if it does not exist.[/dim]",
            title="  Step 5 of 5 · Output folder  ",
            box=box.ROUNDED, border_style="cyan", padding=(0, 2)))
    path = _q(questionary.path,
              "  Output directory (created if it does not exist):",
              default=str(Path.home() / "molscreener_output"),
              only_directories=True)
    return str(Path(path).expanduser().resolve())


def tui_show_summary(cfg: RunConfig) -> bool:
    if console:
        c = cfg.criteria
        t = Table(box=box.SIMPLE_HEAVY, show_header=False, border_style="dim",
                  pad_edge=False, show_edge=False)
        t.add_column("Key",   style="cyan",  no_wrap=True, min_width=24)
        t.add_column("Value", style="white", no_wrap=False)

        # ── input ─────────────────────────────────────────────────────
        t.add_row("[dim]── Libraries ───────────────────────────────────────────[/dim]", "")
        for lib in cfg.libraries:
            t.add_row(f"  {lib.label}", lib.path)

        # ── output ────────────────────────────────────────────────────
        t.add_row("[dim]── Output ────────────────────────────────────────────[/dim]", "")
        t.add_row("  Mode",
                  "3D SDF  →  GNINA / Vina / AutoDock"
                  if cfg.mode == OutputMode.gnina
                  else "Boltz-2 YAMLs  →  boltz predict")
        t.add_row("  Folder", cfg.output_dir)
        if cfg.mode == OutputMode.gnina:
            t.add_row("  Protonation pH", str(cfg.ph))
            t.add_row("  Max molecules/SDF", f"{cfg.batch_size:,}")
        else:
            seq = cfg.boltz2.protein_sequence
            t.add_row("  Receptor sequence",
                      (seq[:50] + f"… ({len(seq)} aa)") if len(seq) > 50
                      else f"{seq}  ({len(seq)} aa)")
            t.add_row("  Receptor MSA",
                      cfg.boltz2.msa_file or "[yellow](none — reduced accuracy)[/yellow]")
            t.add_row("  Chain IDs",
                      f"receptor = {cfg.boltz2.protein_id}   "
                      f"ligand = {cfg.boltz2.ligand_id}")

        # ── compute ───────────────────────────────────────────────────
        t.add_row("[dim]── Compute ────────────────────────────────────────────[/dim]", "")
        t.add_row("  CPUs",       f"{cfg.n_cpus} of {TOTAL_CPUS} available")
        t.add_row("  Chunk size", f"{cfg.chunk_size:,} molecules per task")

        # ── filters ───────────────────────────────────────────────────
        t.add_row("[dim]── Filters ────────────────────────────────────────────[/dim]", "")
        t.add_row("  MW (Da)",           f"{c.min_mw} – {c.max_mw}")
        t.add_row("  LogP",              f"{c.min_logp} – {c.max_logp}")
        t.add_row("  HBD / HBA",         f"{c.min_hbd}–{c.max_hbd}  /  "
                                          f"{c.min_hba}–{c.max_hba}")
        t.add_row("  TPSA (Å²)",         f"{c.min_tpsa} – {c.max_tpsa}")
        t.add_row("  Rotatable bonds",   f"{c.min_rot_bonds} – {c.max_rot_bonds}")
        t.add_row("  fCsp³",             f"{c.min_fcsp3} – {c.max_fcsp3}")
        t.add_row("  QED",               f"{c.min_qed} – {c.max_qed}")
        t.add_row("  SA Score",          f"{c.min_sa} – {c.max_sa}")
        t.add_row("  Heavy atoms",       f"{c.min_heavy} – {c.max_heavy}")
        t.add_row("  Formal charge",     f"{c.min_charge} – {c.max_charge}")
        t.add_row("  Rings / Aromatic",  f"{c.min_rings}–{c.max_rings}  /  "
                                          f"{c.min_arom}–{c.max_arom}")
        t.add_row("  Stereocentres",     f"{c.min_stereo} – {c.max_stereo}")
        t.add_row("  PAINS",
                  "[green]remove[/green]" if c.remove_pains else "[dim]keep[/dim]")
        t.add_row("  Brenk alerts",
                  "[green]remove[/green]" if c.remove_brenk else "[dim]keep[/dim]")

        console.print(Panel(t, title="  Ready to run — please review  ",
                            box=box.ROUNDED, border_style="cyan", padding=(0, 1)))
    else:
        print("\n──── Run Summary ────")
        for lib in cfg.libraries:
            print(f"  Library : {lib.label}  →  {lib.path}")
        print(f"  Mode    : {cfg.mode.value}    CPUs: {cfg.n_cpus}    "
              f"Chunk: {cfg.chunk_size}")
        print(f"  Output  : {cfg.output_dir}")
        c = cfg.criteria
        print(f"  MW {c.min_mw}–{c.max_mw}  LogP {c.min_logp}–{c.max_logp}  "
              f"QED {c.min_qed}–{c.max_qed}")

    return _q(questionary.confirm, "Everything looks good — start screening?",
              default=True)


def tui_save_config(cfg: RunConfig) -> None:
    if console:
        console.print(
            "\n  [dim]Tip: saving to JSON lets you re-run identically with\n"
            "  [cyan]python MolScreener.py --config molscreener_config.json[/cyan]"
            " — useful for SLURM scripts or repeating a campaign.[/dim]")
    if not _q(questionary.confirm,
              "Save this configuration to JSON for future runs?", default=False):
        return
    path = (_q(questionary.text, "  Config filename:",
               default="molscreener_config.json") or "molscreener_config.json")
    data = {
        "libraries":  [{"label": l.label, "path": l.path} for l in cfg.libraries],
        "criteria":   asdict(cfg.criteria),
        "mode":       cfg.mode.value,
        "output_dir": cfg.output_dir,
        "n_cpus":     cfg.n_cpus,
        "chunk_size": cfg.chunk_size,
        "ph":         cfg.ph,
        "batch_size": cfg.batch_size,
        "boltz2":     asdict(cfg.boltz2),
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    _log("ok", f"Config saved → {path}")


def run_tui() -> None:
    if not _Q:
        sys.exit("questionary not installed — pip install questionary")
    if console:
        console.print(Panel(
            "[bold white]MolScreener[/bold white]"
            "  [dim]— compound library preparation for virtual screening[/dim]\n\n"
            "  This tool takes raw vendor libraries and produces clean,\n"
            "  named, filtered, docking-ready inputs in three steps:\n\n"
            "  [bold cyan]1[/bold cyan]  Read your library  "
            "[dim](ZIP · SDF · SDF.gz · SMILES · CSV · folder)[/dim]\n"
            "  [bold cyan]2[/bold cyan]  Filter by drug-like properties  "
            "[dim](MW, LogP, QED, PAINS, Brenk…)[/dim]\n"
            "  [bold cyan]3[/bold cyan]  Write inputs  "
            "[dim]— 3D SDF ready for GNINA / Vina  "
            "or  YAML ready for Boltz-2[/dim]\n\n"
            "  [dim]Every compound keeps its catalogue ID throughout.\n"
            "  Use Ctrl-C at any time to stop — partial output is always saved.[/dim]",
            title="  Welcome  ",
            box=box.DOUBLE_EDGE, border_style="cyan", padding=(1, 2)))
    else:
        print("\n  MolScreener — compound library preparation for virtual screening")
        print("  Library → filter → 3D SDF (GNINA) or YAML (Boltz-2)\n")

    libs = tui_collect_libraries()
    crit = tui_collect_criteria()
    mode, ph, batch, bcfg = tui_collect_mode()
    n_cpus, chunk = tui_collect_compute()
    outdir = tui_collect_output()

    cfg = RunConfig(libraries=libs, criteria=crit, mode=mode,
                    output_dir=outdir, n_cpus=n_cpus, chunk_size=chunk,
                    ph=ph, batch_size=batch, boltz2=bcfg or Boltz2Config())

    if not tui_show_summary(cfg):
        _log("warn", "Aborted by user."); return
    tui_save_config(cfg)

    stop_flag = [False]
    try:
        run_pipeline(cfg, stop_flag)
    except KeyboardInterrupt:
        stop_flag[0] = True
        _erase_progress()
        _log("warn", "Interrupted — partial output saved.")


# ══════════════════════════════════════════════════════════════════════
# §11  CLI
# ══════════════════════════════════════════════════════════════════════

def run_cli(args) -> None:
    if args.config:
        with open(args.config) as f:
            data = json.load(f)
        libs = [LibraryEntry(**l) for l in data.get("libraries", [])]
        crit_d = data.get("criteria", {})
        crit = ScreeningCriteria(**crit_d) if crit_d else ScreeningCriteria()
        mode = OutputMode(data.get("mode", "gnina"))
        bcfg = Boltz2Config(**data.get("boltz2", {}))
        cfg = RunConfig(
            libraries=libs, criteria=crit, mode=mode,
            output_dir=data.get("output_dir", ""),
            n_cpus=int(data.get("n_cpus", max(1, TOTAL_CPUS - 1))),
            chunk_size=int(data.get("chunk_size", 2000)),
            ph=float(data.get("ph", 7.0)),
            batch_size=int(data.get("batch_size", 100_000)),
            boltz2=bcfg)
    else:
        libs = []
        for entry in (args.lib or []):
            path, label = entry[0], entry[1]
            if not os.path.exists(path):
                sys.exit(f"Library path not found: {path}")
            if not resolve_shards(path):
                sys.exit(f"No .sdf/.smi/.csv content found in: {path}")
            libs.append(LibraryEntry(label=label, path=path))
        if not libs:
            sys.exit("No libraries specified. Use --lib PATH LABEL or --config.")
        cfg = RunConfig(
            libraries=libs,
            mode=OutputMode(args.mode),
            output_dir=args.output or "molscreener_output",
            n_cpus=max(1, min(args.cpus or max(1, TOTAL_CPUS - 1), TOTAL_CPUS)),
            chunk_size=max(100, args.chunk),
            ph=args.ph, batch_size=args.batch)
        if cfg.mode == OutputMode.boltz2:
            if not args.sequence:
                sys.exit("--sequence is required for boltz2 mode")
            msa = (str(Path(args.msa).expanduser().resolve())
                   if args.msa else "")
            cfg.boltz2 = Boltz2Config(
                protein_sequence=args.sequence, msa_file=msa,
                protein_id=args.protein_id or "A",
                ligand_id=args.ligand_id or "B")

    if not cfg.output_dir:
        sys.exit("output_dir is required.")
    for lib in cfg.libraries:
        if not os.path.exists(lib.path):
            sys.exit(f"Library not found: {lib.path}")

    stop_flag = [False]
    try:
        run_pipeline(cfg, stop_flag)
    except KeyboardInterrupt:
        stop_flag[0] = True
        _erase_progress()
        _log("warn", "Interrupted — partial output saved.")


# ══════════════════════════════════════════════════════════════════════
# §12  ENTRY POINT
# ══════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="MolScreener",
        description="Compound library preparation for virtual screening — "
                    "ZIP/SDF/SMILES/CSV  →  3D SDF (GNINA) or YAML (Boltz-2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
────────
  # Interactive TUI (recommended for first use):
  python MolScreener.py

  # Re-run a saved configuration (great for SLURM scripts):
  python MolScreener.py --config molscreener_config.json

  # Mixed libraries, GNINA docking prep, 32 CPUs:
  python MolScreener.py \\
    --lib /data/Enamine.sdf.gz  Enamine \\
    --lib /data/MolPort.zip     MolPort \\
    --lib /data/custom.smi      Custom  \\
    --mode gnina --cpus 32 --chunk 2000 --ph 7.4 \\
    --output /scratch/gnina_ready

  # Boltz-2 affinity prep, directory of SDFs, 16 CPUs:
  python MolScreener.py \\
    --lib /data/library/        MyLib   \\
    --mode boltz2 --cpus 16    \\
    --sequence MTEYKLVVVGAGG...         \\
    --msa /data/GHRHR_UNPAIRED_MSA.a3m \\
    --output /scratch/boltz2_ready
""")
    parser.add_argument("--config", metavar="JSON",
        help="Path to a saved JSON config (skips all other flags)")
    parser.add_argument("--lib", metavar=("PATH", "LABEL"), nargs=2, action="append",
        help="Library path (ZIP/SDF/SMILES/CSV/dir) + label (repeatable)")
    parser.add_argument("--zip", metavar=("PATH", "LABEL"), nargs=2, action="append",
        dest="lib", help="Alias for --lib (back-compat)")
    parser.add_argument("--mode", choices=["gnina", "boltz2"], default="gnina")
    parser.add_argument("--output", metavar="DIR", help="Output directory")
    parser.add_argument("--cpus", type=int, help=f"CPUs to use (max {TOTAL_CPUS})")
    parser.add_argument("--chunk", type=int, default=500,
        help="Molecules per parallel task (default 500)")
    parser.add_argument("--ph", type=float, default=7.0,
        help="Protonation pH for GNINA mode (default 7.0)")
    parser.add_argument("--batch", type=int, default=100_000,
        help="Molecules per output SDF file (GNINA, default 100000)")
    parser.add_argument("--sequence", metavar="SEQ",
        help="Protein FASTA sequence — required for boltz2 mode")
    parser.add_argument("--msa", metavar="FILE", help="MSA .a3m file for Boltz-2")
    parser.add_argument("--protein-id", default="A")
    parser.add_argument("--ligand-id", default="B")

    args = parser.parse_args()
    if len(sys.argv) == 1:
        run_tui()
    else:
        run_cli(args)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    # On macOS (Python ≥ 3.8) the default start method is 'spawn', which
    # requires this guard to be present — otherwise worker processes
    # re-import the script and silently fall back to a single CPU.
    # Setting it explicitly here ensures correct behaviour on macOS, Linux,
    # and HPC clusters regardless of the system default.
    multiprocessing.set_start_method("spawn", force=True)
    main()
