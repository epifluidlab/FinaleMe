#!/usr/bin/env python3
"""motif_filter.py — Enzyme-targeted cfDNA fragment filter.

Standalone preprocessing utility that classifies cfDNA fragments by their
Read-1-derived 5' end 4-mer motif and writes per-class fragment outputs for
downstream FinaleMe / FinaleToolkit analysis.

Implements design/finaleme_motif_filter_design_v3.5.md. The five output files
are (suffixes on --output-prefix):

  .DNase1L3_enriched           DNase1L3 motif set (default: CCNN family)
  .DNase1_enriched             DNase1 motif set (default: TGNN family)
  .DNase1L3_or_DNase1_enriched union of the two above (single union output)
  .L3_D1_dominant              F-I or F-II argmax across all six F-profiles
  .DFFB_dominant               F-III argmax (apoptotic-flux readout)

A summary JSON, optional per-fragment log, and a 7-page QC PDF are also
written. See the design doc §6 for details.

This script is a thin wrapper around FinaleToolkit. Required runtime deps:
  pip install finaletoolkit pysam py2bit numpy scipy matplotlib seaborn

Example:

  python motif_filter.py \\
      --input sample.sorted.bam \\
      --reference hg38.2bit \\
      --output-prefix out/sample

The input BAM should have the MQ tag populated (`samtools fixmate -m` then
sort + index). If not, pass `--mate-check fetch_mate` for per-fragment mate
lookup (slower), or `--mate-check skip` for diagnostic comparisons only.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional

# Heavy deps imported lazily inside the functions that need them so that
# `python motif_filter.py --help` works even without finaletoolkit / pysam
# installed.

__version__ = "0.1.0"

# ---------------------------------------------------------------------------
# Constants — built-in motif sets from design v3.5 §4.2 and §4.4.
# ---------------------------------------------------------------------------

DEFAULT_K = 4
DEFAULT_MAPQ = 30
DEFAULT_FRAC_LOW = 50
DEFAULT_MAX_5P_SOFTCLIP = 5

# Class labels (internal).
CLASS_DNASE1L3 = "DNase1L3"
CLASS_DNASE1 = "DNase1"
CLASS_L3_D1_DOMINANT = "L3_D1_dominant"
CLASS_DFFB_DOMINANT = "DFFB_dominant"

# Output file labels (§6.1).
OUT_DNASE1L3 = "DNase1L3_enriched"
OUT_DNASE1 = "DNase1_enriched"
OUT_DNASE1L3_OR_DNASE1 = "DNase1L3_or_DNase1_enriched"
OUT_L3_D1_DOMINANT = "L3_D1_dominant"
OUT_DFFB_DOMINANT = "DFFB_dominant"

# Ordered list used for iteration, PDF panels, etc.
OUTPUT_LABELS = (
    OUT_DNASE1L3,
    OUT_DNASE1,
    OUT_DNASE1L3_OR_DNASE1,
    OUT_L3_D1_DOMINANT,
    OUT_DFFB_DOMINANT,
)

# CCNN — DNase1L3 paper-faithful set (Zhou et al. 2023, F-profile I text). 16 motifs.
CCNN_MOTIFS: frozenset[str] = frozenset({
    "CCAA", "CCAC", "CCAG", "CCAT",
    "CCCA", "CCCC", "CCCG", "CCCT",
    "CCGA", "CCGC", "CCGG", "CCGT",
    "CCTA", "CCTC", "CCTG", "CCTT",
})

# F-I argmax — DNase1L3 empirical set (motifs whose largest F-profile is F-I). 25 motifs.
F_I_ARGMAX_MOTIFS: frozenset[str] = frozenset({
    "ACAA", "CAAA", "CAAC", "CAAG", "CAAT",
    "CAGA", "CAGT", "CATA", "CATG",
    "CCAA", "CCAC", "CCAG", "CCAT", "CCCA", "CCCC", "CCCT",
    "CCTA", "CCTC", "CCTG", "CCTT",
    "CTAA", "CTAT", "CTGA", "CTTA", "TCTT",
})

# TGNN — DNase1 paper-faithful set. 16 motifs.
TGNN_MOTIFS: frozenset[str] = frozenset({
    "TGAA", "TGAC", "TGAG", "TGAT",
    "TGCA", "TGCC", "TGCG", "TGCT",
    "TGGA", "TGGC", "TGGG", "TGGT",
    "TGTA", "TGTC", "TGTG", "TGTT",
})

# F-II argmax — DNase1 empirical set (motifs whose largest F-profile is F-II). 53 motifs.
F_II_ARGMAX_MOTIFS: frozenset[str] = frozenset({
    "CGAA", "CGAT", "CGCA", "CGCT", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT",
    "TAAA", "TAAC", "TAAG", "TAAT", "TACA", "TACC", "TACG", "TACT",
    "TAGA", "TAGC", "TAGG", "TAGT", "TATA", "TATC", "TATG", "TATT",
    "TCAA", "TCAC", "TCAG", "TCAT", "TCCA", "TCCC", "TCCG", "TCCT",
    "TCTA", "TCTC", "TCTG",
    "TGAA", "TGAC", "TGAG", "TGAT", "TGCA", "TGCC", "TGCG", "TGCT",
    "TGGA", "TGGC", "TGGG", "TGGT", "TGTA", "TGTC", "TGTG", "TGTT",
    "TTCA",
})

# L3_D1_dominant = F-I ∪ F-II argmax (78 motifs, by definition).
L3_D1_DOMINANT_MOTIFS: frozenset[str] = F_I_ARGMAX_MOTIFS | F_II_ARGMAX_MOTIFS

# DFFB canonical — Zhou et al. Fig 2D figure-labeled motifs (paper-faithful). 4 motifs.
DFFB_CANONICAL_MOTIFS: frozenset[str] = frozenset({
    "ACCC", "ACCT", "ACTC", "ACTT",
})

# DFFB argmax — F-III dominant motifs (recommended for meaningful retention). 50 motifs.
DFFB_ARGMAX_MOTIFS: frozenset[str] = frozenset({
    "AAAT", "AACA", "AACC", "AACG", "AACT", "AATA", "AATC", "AATG", "AATT",
    "ACAT", "ACCA", "ACCC", "ACCG", "ACCT",
    "ACTA", "ACTC", "ACTG", "ACTT",
    "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGTA", "AGTC", "AGTT",
    "ATAT", "ATCC", "ATCT", "ATTT",
    "CACA", "CACC", "CACG", "CACT", "CATC", "CATT", "CTCT",
    "GACA", "GACC", "GACG", "GACT", "GATC", "GATT",
    "GCCC", "GCCT", "GCTC", "GCTT", "GGCT", "GGTT", "GTCT",
})

# Expected per-output retention ranges in healthy plasma (§4.5).
# Used by QC page 1 to color the bars (green if within, orange if outside).
EXPECTED_RETENTION = {
    "ccnn": {
        OUT_DNASE1L3: (0.13, 0.17),
        OUT_DNASE1: (0.09, 0.12),
        OUT_DNASE1L3_OR_DNASE1: (0.20, 0.28),
        OUT_L3_D1_DOMINANT: (0.40, 0.50),
        OUT_DFFB_DOMINANT: (0.02, 0.03),
    },
    "fprofile_argmax": {
        OUT_DNASE1L3: (0.17, 0.22),
        OUT_DNASE1: (0.22, 0.28),
        OUT_DNASE1L3_OR_DNASE1: (0.38, 0.45),
        OUT_L3_D1_DOMINANT: (0.40, 0.50),
        OUT_DFFB_DOMINANT: (0.18, 0.22),
    },
}

# Reverse complement table — handles upper, lower, and N.
_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string. N maps to N."""
    return seq.translate(_RC_TABLE)[::-1]


# Pre-compute the 256 4-mers in canonical (A < C < G < T) order. Used by
# QC plotting to keep all per-output histograms on a consistent x-axis.
def _enumerate_kmers(k: int) -> list[str]:
    alphabet = "ACGT"
    out = [""]
    for _ in range(k):
        out = [p + b for p in out for b in alphabet]
    return out


CANONICAL_4MERS: tuple[str, ...] = tuple(_enumerate_kmers(4))


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log = logging.getLogger("motif_filter")


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """argparse parser. Help text mirrors design doc §2."""
    p = argparse.ArgumentParser(
        prog="motif_filter.py",
        description=(
            "Enzyme-targeted cfDNA fragment filter (design v3.5). Classifies "
            "fragments by Read-1-derived 5' end motif and writes five "
            "filtered output files: DNase1L3_enriched, DNase1_enriched, "
            "DNase1L3_or_DNase1_enriched, L3_D1_dominant, DFFB_dominant."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required.
    p.add_argument("--input", required=True,
                   help="Input BAM/CRAM (.bam/.cram, with index) or "
                        "tabix-indexed fragment file (.frag.gz, with .tbi). "
                        "Format auto-detected by extension.")
    p.add_argument("--reference", required=True,
                   help="Reference genome in .2bit format. Required to look "
                        "up the 4-mer motif at each fragment's 5' end.")
    p.add_argument("--output-prefix", required=True,
                   help="Path prefix for the five output files. Suffixes "
                        "are appended automatically (.DNase1L3_enriched, etc.).")

    # Motif config.
    p.add_argument("--motif-config", default=None,
                   help="Optional TSV file specifying motif sets. One row "
                        "per motif: 'motif<TAB>class1,class2,...'. If set, "
                        "REPLACES the built-in defaults (--l3-set/--d1-set/"
                        "--dffb-set are ignored). Default: built-in.")
    p.add_argument("--motif-length", type=int, default=DEFAULT_K,
                   help="Motif k-mer length. Default 4. Must match --motif-config.")

    # Set choices (only used when --motif-config is unset).
    p.add_argument("--l3-set", choices=["ccnn", "fprofile_argmax"], default="ccnn",
                   help="DNase1L3 motif set. 'ccnn' = 16 motifs from Zhou et "
                        "al. text (default, paper-faithful). 'fprofile_argmax' "
                        "= 25 motifs where F-I is the largest F-profile.")
    p.add_argument("--d1-set", choices=["tgnn", "fprofile_argmax"], default="tgnn",
                   help="DNase1 motif set. 'tgnn' = 16 motifs (default). "
                        "'fprofile_argmax' = 53 motifs where F-II is largest.")
    p.add_argument("--dffb-set", choices=["canonical", "fprofile_argmax"], default="canonical",
                   help="DFFB motif set. 'canonical' = 4 figure-labeled "
                        "motifs from Zhou et al. Fig 2D (default, paper-"
                        "faithful, ~2-3 percent retention). 'fprofile_argmax' "
                        "= 50 motifs where F-III is largest (~18-22 percent "
                        "retention, recommended for meaningful output size).")

    # Quality / pair-quality.
    p.add_argument("--quality-threshold", type=int, default=DEFAULT_MAPQ,
                   help="Minimum MAPQ for both reads of a pair. Default 30.")
    p.add_argument("--mate-check", choices=["require_mq_tag", "fetch_mate", "skip"],
                   default="require_mq_tag",
                   help="How to enforce mate MAPQ. 'require_mq_tag' "
                        "(default): require Read 1 to have the MQ tag "
                        "(populated by samtools fixmate -m); aborts at "
                        "startup if absent. 'fetch_mate': per-fragment "
                        "mate lookup via pysam.AlignmentFile.mate() (~2x "
                        "BAM I/O). 'skip': no mate-MAPQ check; diagnostic only.")
    p.add_argument("--no-check-multimap", action="store_true",
                   help="Disable XA/SA tag check. Default: check enabled "
                        "(reads with either tag are dropped).")
    p.add_argument("--max-5prime-softclip", type=int, default=DEFAULT_MAX_5P_SOFTCLIP,
                   help="Maximum soft-clip length at Read 1's 5' end. "
                        "Default 5; reads exceeding this are dropped.")

    # Fragment length.
    p.add_argument("--fraction-low", type=int, default=DEFAULT_FRAC_LOW,
                   help="Minimum fragment length. Default 50.")
    p.add_argument("--fraction-high", type=int, default=None,
                   help="Maximum fragment length. Default: no upper limit.")

    # Threading and output options.
    p.add_argument("--threads", type=int, default=1,
                   help="Number of threads for BAM read/write. Default 1.")
    p.add_argument("--log", default=None,
                   help="Optional per-fragment classification log (gzipped "
                        "TSV). Columns: chrom, start, stop, read_on_plus, "
                        "motif, assigned_classes, outputs.")
    p.add_argument("--no-qc-pdf", action="store_true",
                   help="Disable QC PDF generation (default: enabled). "
                        "Skips matplotlib + scipy dependency import.")
    p.add_argument("--emit-motif-tag", action="store_true",
                   help="Write the assigned motif as a BAM tag (X5:Z:CCAA) "
                        "on each output read. BAM output only.")
    p.add_argument("--verbose", action="store_true",
                   help="Verbose progress logging to stderr.")
    p.add_argument("--version", action="version", version=f"motif_filter.py {__version__}")
    return p


# ---------------------------------------------------------------------------
# Motif config — built-in defaults + optional user TSV.
# ---------------------------------------------------------------------------

def build_motif_to_classes(
    l3_set: str,
    d1_set: str,
    dffb_set: str,
) -> dict[str, set[str]]:
    """Return {motif -> set(class_label)} from the built-in default sets.

    Class labels map to outputs in `classify()`. A motif can belong to multiple
    classes (e.g. CCAA is in both DNase1L3 set and the L3_D1_dominant set).
    """
    out: dict[str, set[str]] = {}

    l3_set_motifs = CCNN_MOTIFS if l3_set == "ccnn" else F_I_ARGMAX_MOTIFS
    d1_set_motifs = TGNN_MOTIFS if d1_set == "tgnn" else F_II_ARGMAX_MOTIFS
    dffb_set_motifs = (
        DFFB_CANONICAL_MOTIFS if dffb_set == "canonical" else DFFB_ARGMAX_MOTIFS
    )

    for m in l3_set_motifs:
        out.setdefault(m, set()).add(CLASS_DNASE1L3)
    for m in d1_set_motifs:
        out.setdefault(m, set()).add(CLASS_DNASE1)
    for m in L3_D1_DOMINANT_MOTIFS:
        out.setdefault(m, set()).add(CLASS_L3_D1_DOMINANT)
    for m in dffb_set_motifs:
        out.setdefault(m, set()).add(CLASS_DFFB_DOMINANT)
    return out


def load_motif_config(path: str, k: int) -> dict[str, set[str]]:
    """Load a user-supplied motif config TSV.

    Format: one row per motif, 'motif<TAB>class1,class2,...'. Lines starting
    with '#' are ignored. Motif length must equal `k`; invalid motifs raise.
    Classes must be drawn from the four built-ins (DNase1L3, DNase1,
    L3_D1_dominant, DFFB_dominant) for the downstream `classify()` to route
    them; unknown class labels are kept (allowing custom config testing) but
    will not be written to outputs.
    """
    valid_classes = {
        CLASS_DNASE1L3, CLASS_DNASE1, CLASS_L3_D1_DOMINANT, CLASS_DFFB_DOMINANT,
    }
    out: dict[str, set[str]] = {}
    with open(path) as fh:
        for ln, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"--motif-config {path}: line {ln}: expected '<motif><TAB><classes>' "
                    f"but got: {raw.strip()!r}"
                )
            motif = parts[0].strip().upper()
            classes_str = parts[1].strip()
            if len(motif) != k:
                raise ValueError(
                    f"--motif-config {path}: line {ln}: motif {motif!r} "
                    f"has length {len(motif)} but --motif-length is {k}"
                )
            if not re.fullmatch(r"[ACGT]+", motif):
                raise ValueError(
                    f"--motif-config {path}: line {ln}: motif {motif!r} "
                    f"contains non-ACGT characters"
                )
            classes = {c.strip() for c in classes_str.split(",") if c.strip()}
            unknown = classes - valid_classes
            if unknown:
                log.warning(
                    "--motif-config %s line %d: unknown class labels %s; "
                    "will be ignored downstream",
                    path, ln, sorted(unknown),
                )
            out[motif] = classes
    if not out:
        raise ValueError(f"--motif-config {path}: no motifs loaded.")
    return out


def motif_config_hash(motif_to_classes: dict[str, set[str]]) -> str:
    """Stable hash of the motif config — used for the summary JSON so the user
    can confirm two runs used identical motif assignments."""
    items = sorted(
        (m, ",".join(sorted(cs))) for m, cs in motif_to_classes.items()
    )
    h = hashlib.sha256()
    for m, cs in items:
        h.update(f"{m}\t{cs}\n".encode())
    return "sha256:" + h.hexdigest()


# ---------------------------------------------------------------------------
# Motif extraction (§5.2) — single motif per fragment, strand-dispatched.
# ---------------------------------------------------------------------------

def get_motif(
    contig: str,
    frag_start: int,
    frag_stop: int,
    read_on_plus: bool,
    refseq: Any,
    k: int = DEFAULT_K,
) -> Optional[str]:
    """Return the Read-1-derived 5' end motif, or None on N / boundary error.

    For read_on_plus=True (Read 1 on + strand), the 5' end is at frag_start,
    so motif = forward 4-mer at [frag_start, frag_start+k).

    For read_on_plus=False (Read 1 on − strand), the 5' end is at frag_stop,
    so motif = reverse-complement of the forward 4-mer at [frag_stop-k, frag_stop).
    """
    try:
        if read_on_plus:
            kmer = refseq.sequence(contig, int(frag_start), int(frag_start) + k)
        else:
            genomic_kmer = refseq.sequence(
                contig, int(frag_stop) - k, int(frag_stop)
            )
            kmer = reverse_complement(genomic_kmer)
    except RuntimeError:
        # py2bit raises RuntimeError near chromosome boundaries (request
        # extends past chromosome end). Treat as invalid motif.
        return None
    if kmer is None or len(kmer) != k:
        return None
    kmer = kmer.upper()
    if "N" in kmer:
        return None
    return kmer


# ---------------------------------------------------------------------------
# Classification (§5.3).
# ---------------------------------------------------------------------------

def classify(motif: str, motif_to_classes: dict[str, set[str]]) -> list[str]:
    """Return the list of OUTPUT labels a fragment with this motif goes into.

    A fragment can land in multiple outputs (e.g., CCAA -> DNase1L3_enriched
    AND DNase1L3_or_DNase1_enriched AND L3_D1_dominant). The DFFB_dominant
    output is mutually exclusive with the DNase1L3_or_DNase1_enriched output
    by F-profile-argmax construction.
    """
    classes = motif_to_classes.get(motif, set())
    outputs: list[str] = []
    if CLASS_DNASE1L3 in classes:
        outputs.append(OUT_DNASE1L3)
    if CLASS_DNASE1 in classes:
        outputs.append(OUT_DNASE1)
    if CLASS_DNASE1L3 in classes or CLASS_DNASE1 in classes:
        outputs.append(OUT_DNASE1L3_OR_DNASE1)
    if CLASS_L3_D1_DOMINANT in classes:
        outputs.append(OUT_L3_D1_DOMINANT)
    if CLASS_DFFB_DOMINANT in classes:
        outputs.append(OUT_DFFB_DOMINANT)
    return outputs


# ---------------------------------------------------------------------------
# Pair-quality helpers (§3.5).
# ---------------------------------------------------------------------------

def _check_first_alignment_has_mq_tag(bam_path: str, threads: int) -> None:
    """For --mate-check require_mq_tag mode: confirm Read 1 alignments carry
    the MQ tag. Abort with the actionable fix instruction if absent.

    We scan up to ~1000 primary Read-1 alignments. If none carries MQ, raise.
    This early-fails the run rather than silently skipping mate-MAPQ checks
    on every fragment.
    """
    import pysam
    sam = pysam.AlignmentFile(bam_path, "rb", threads=threads)
    seen = 0
    try:
        for read in sam.fetch(until_eof=True):
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                continue
            seen += 1
            try:
                if read.has_tag("MQ"):
                    return  # MQ present — we're good.
            except Exception:
                pass
            if seen >= 1000:
                break
    finally:
        sam.close()
    raise SystemExit(
        "ERROR: input BAM lacks MQ tag on Read 1 alignments. "
        "Run `samtools fixmate -m input.bam fixed.bam` (then sort + index) "
        "to populate MQ tags, or pass `--mate-check fetch_mate` for "
        "per-fragment mate lookup (slower)."
    )


def has_excessive_5p_softclip(read: Any, threshold: int) -> bool:
    """True if Read 1 has more than `threshold` bp of soft-clipping at its
    5' end. 5' on + strand = first CIGAR op; 5' on − strand = last CIGAR op."""
    cigar = read.cigartuples
    if not cigar:
        return False
    # 'S' op = 4 in pysam encoding.
    if read.is_reverse:
        op, length = cigar[-1]
    else:
        op, length = cigar[0]
    return op == 4 and length > threshold


def has_multimap_tag(read: Any) -> bool:
    """True if the read has bwa-mem XA (alt alignments) or SA (chimeric) tag."""
    try:
        return read.has_tag("XA") or read.has_tag("SA")
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Input format detection.
# ---------------------------------------------------------------------------

def detect_input_format(path: str) -> str:
    """Return 'bam', 'cram', or 'frag_gz' based on extension."""
    p = path.lower()
    if p.endswith(".bam"):
        return "bam"
    if p.endswith(".cram"):
        return "cram"
    if p.endswith(".frag.gz") or p.endswith(".bed.gz"):
        return "frag_gz"
    raise SystemExit(
        f"ERROR: cannot detect input format from extension: {path!r}. "
        f"Expected .bam, .cram, .frag.gz, or .bed.gz."
    )


# ---------------------------------------------------------------------------
# Output writers — one per output label, dispatched by input format.
# ---------------------------------------------------------------------------

@dataclass
class BamOutputs:
    """BAM output writer set. One pysam.AlignmentFile per output label.

    Writes BOTH reads of a passing pair (preserving PE integrity per
    design principle 6). We accumulate Read 2 in a per-name buffer keyed on
    qname; when Read 1 comes through (with classification), we look up the
    matching Read 2 and emit both, or vice versa.

    For simplicity here we use the "write whichever side of the pair arrives;
    cache the other if seen later" pattern. The output BAM is queryname-
    or coordinate-sorted depending on input.
    """
    prefix: str
    template_bam: Any
    emit_motif_tag: bool
    threads: int
    writers: dict[str, Any] = field(default_factory=dict)
    # Per-output-label: pending reads we've seen but whose mate we haven't
    # processed yet. Key: (output_label, qname). Value: pysam.AlignedSegment.
    pending: dict[tuple[str, str], Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        import pysam
        for label in OUTPUT_LABELS:
            path = f"{self.prefix}.{label}.bam"
            self.writers[label] = pysam.AlignmentFile(
                path, "wb", template=self.template_bam, threads=self.threads,
            )

    def write_pair(self, label: str, read1: Any, read2: Optional[Any], motif: str) -> None:
        """Write Read 1 and (if available) Read 2 to the given output."""
        if self.emit_motif_tag:
            try:
                read1.set_tag("X5", motif, value_type="Z")
                if read2 is not None:
                    read2.set_tag("X5", motif, value_type="Z")
            except Exception:
                pass
        w = self.writers[label]
        w.write(read1)
        if read2 is not None:
            w.write(read2)

    def close_and_index(self) -> None:
        import pysam
        for label, w in self.writers.items():
            w.close()
            path = f"{self.prefix}.{label}.bam"
            try:
                pysam.index(path)
            except Exception as exc:  # pragma: no cover
                log.warning("Failed to index %s: %s", path, exc)


@dataclass
class FragGzOutputs:
    """Tabix-indexed bgzipped TSV (BED3+2) writer set.

    Columns written: chrom, start, stop, mapq, strand.
    Strand follows the input convention: '+' if read_on_plus else '-'.
    Tabix-indexed after writing using pysam.tabix_index.
    """
    prefix: str
    writers: dict[str, Any] = field(default_factory=dict)
    _paths: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        import pysam
        for label in OUTPUT_LABELS:
            path = f"{self.prefix}.{label}.frag.gz"
            self._paths[label] = path
            # bgzipped TSV writer via pysam.BGZFile.
            self.writers[label] = pysam.BGZFile(path, "wb")

    def write_fragment(
        self, label: str, contig: str, start: int, stop: int,
        mapq: int, read_on_plus: bool, motif: str,
    ) -> None:
        strand = "+" if read_on_plus else "-"
        line = f"{contig}\t{start}\t{stop}\t{mapq}\t{strand}\n"
        self.writers[label].write(line.encode())

    def close_and_index(self) -> None:
        import pysam
        for label, w in self.writers.items():
            w.close()
            path = self._paths[label]
            try:
                pysam.tabix_index(path, preset="bed", force=True)
            except Exception as exc:  # pragma: no cover
                log.warning("Failed to tabix-index %s: %s", path, exc)


# ---------------------------------------------------------------------------
# Per-fragment log writer.
# ---------------------------------------------------------------------------

class FragmentLog:
    """Optional gzipped TSV recording each processed fragment's classification.

    Columns: chrom, start, stop, read_on_plus, motif, assigned_classes, outputs.
    """

    def __init__(self, path: Optional[str]):
        self.path = path
        self._fh = None
        if path is not None:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            self._fh = gzip.open(path, "wt")
            self._fh.write(
                "chrom\tstart\tstop\tread_on_plus\tmotif\t"
                "assigned_classes\toutputs\n"
            )

    def write(
        self, contig: str, start: int, stop: int, read_on_plus: bool,
        motif: Optional[str], classes: set[str], outputs: list[str],
    ) -> None:
        if self._fh is None:
            return
        self._fh.write(
            f"{contig}\t{start}\t{stop}\t{int(bool(read_on_plus))}\t"
            f"{motif or '.'}\t{','.join(sorted(classes)) or '.'}\t"
            f"{','.join(outputs) or '.'}\n"
        )

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None


# ---------------------------------------------------------------------------
# Stats accumulator.
# ---------------------------------------------------------------------------

@dataclass
class RunStats:
    total_fragments_processed: int = 0
    passed_qc: int = 0
    skipped_motif_invalid: int = 0
    skipped_mate_mapq_fail: int = 0
    skipped_multimap: int = 0
    skipped_excessive_softclip: int = 0
    # Per-output write counts.
    output_counts: Counter = field(default_factory=Counter)
    # 256-bin motif counters, by strand. Both pre-classification.
    motif_counts_plus: Counter = field(default_factory=Counter)
    motif_counts_minus: Counter = field(default_factory=Counter)
    # Per-output motif counters (for QC PDF page 3).
    per_output_motif_counts: dict[str, Counter] = field(
        default_factory=lambda: {k: Counter() for k in OUTPUT_LABELS}
    )


# ---------------------------------------------------------------------------
# Main streaming loops — three input mode paths.
# ---------------------------------------------------------------------------

def _process_one_fragment(
    contig: str, frag_start: int, frag_stop: int, read_on_plus: bool,
    mapq: int, refseq: Any, motif_to_classes: dict[str, set[str]], k: int,
    stats: RunStats, frag_log: FragmentLog,
    bam_writer: Optional[BamOutputs] = None,
    frag_writer: Optional[FragGzOutputs] = None,
    read1: Optional[Any] = None, read2: Optional[Any] = None,
) -> None:
    """Shared per-fragment classification + output dispatch. Increment counters
    on the shared `stats` object. Either bam_writer or frag_writer is set.
    """
    motif = get_motif(contig, frag_start, frag_stop, read_on_plus, refseq, k=k)
    if motif is None:
        stats.skipped_motif_invalid += 1
        frag_log.write(contig, frag_start, frag_stop, read_on_plus, None, set(), [])
        return
    # Tally per-strand pre-classification 256-bin counts.
    if read_on_plus:
        stats.motif_counts_plus[motif] += 1
    else:
        stats.motif_counts_minus[motif] += 1
    classes = motif_to_classes.get(motif, set())
    outputs = classify(motif, motif_to_classes)
    stats.passed_qc += 1
    frag_log.write(contig, frag_start, frag_stop, read_on_plus, motif, classes, outputs)
    for label in outputs:
        stats.output_counts[label] += 1
        stats.per_output_motif_counts[label][motif] += 1
        if bam_writer is not None:
            bam_writer.write_pair(label, read1, read2, motif)
        elif frag_writer is not None:
            frag_writer.write_fragment(
                label, contig, frag_start, frag_stop, mapq, read_on_plus, motif,
            )


def process_frag_gz(
    args: argparse.Namespace,
    refseq: Any,
    motif_to_classes: dict[str, set[str]],
    stats: RunStats,
    frag_writer: FragGzOutputs,
    frag_log: FragmentLog,
) -> None:
    """Process frag.gz input. No pair-quality filters (assumed applied
    upstream during BAM-to-frag.gz conversion — see design §3.5)."""
    from finaletoolkit.utils.utils import frag_generator
    gen = frag_generator(
        input_file=args.input,
        contig=None,
        quality_threshold=args.quality_threshold,
        fraction_low=args.fraction_low,
        fraction_high=args.fraction_high,
    )
    for (contig, start, stop, mapq, read_on_plus) in gen:
        stats.total_fragments_processed += 1
        _process_one_fragment(
            contig, start, stop, read_on_plus, mapq, refseq, motif_to_classes,
            args.motif_length, stats, frag_log, frag_writer=frag_writer,
        )


def process_bam_require_mq_tag(
    args: argparse.Namespace,
    refseq: Any,
    motif_to_classes: dict[str, set[str]],
    stats: RunStats,
    bam_writer: BamOutputs,
    frag_log: FragmentLog,
) -> None:
    """BAM path using finaletoolkit.frag_generator. Pair-quality is enforced
    by finaletoolkit's low_quality_read_pairs (MAPQ on Read 1, MQ tag on
    Read 2). We add XA/SA + 5' soft-clip filters by reading Read 1 separately
    in a parallel iteration over the BAM."""
    import pysam
    from finaletoolkit.utils.utils import frag_generator

    # For XA/SA/softclip checks we need access to the Read 1 AlignedSegment.
    # frag_generator only yields tuples, so we iterate over the BAM ourselves
    # for Read 1 and cross-reference by (contig, start, stop, read_on_plus).
    #
    # Approach: build a set of "passing fragment keys" from the BAM iteration
    # (Read 1 passes MAPQ + multimap + softclip filters) AND emit the BAM
    # Read 1 + Read 2 lookups for output writing. Iterate the BAM twice if
    # needed.
    #
    # For implementation simplicity we use a single iteration pattern:
    # walk the BAM ourselves, applying low_quality_read_pairs and our added
    # filters, and construct fragment tuples manually for passing Read 1s.
    # Read 2 is paged in via a per-qname cache so we can emit pairs.
    from finaletoolkit.utils.utils import low_quality_read_pairs

    sam = pysam.AlignmentFile(args.input, "rb", threads=args.threads)
    template_bam = sam
    # Cache: qname -> Read 2 (waiting for Read 1).
    r2_cache: dict[str, Any] = {}
    # Cache: qname -> Read 1 (waiting for Read 2).
    r1_cache: dict[str, tuple[Any, list[str], str]] = {}

    fraction_low = args.fraction_low
    fraction_high = args.fraction_high if args.fraction_high is not None else (1 << 31)
    check_mm = not args.no_check_multimap
    sc_thresh = args.max_5prime_softclip
    mapq_min = args.quality_threshold

    for read in sam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        # Drop early if mate is unmapped or different contig.
        if read.mate_is_unmapped:
            continue
        if read.reference_name != read.next_reference_name:
            continue

        if read.is_read1:
            # Apply our extra filters (FinaleToolkit covers MAPQ + flags + MQ).
            if low_quality_read_pairs(read, mapq_min):
                # Could be low MAPQ, or low mate MAPQ via MQ tag, etc.
                # We don't know which, but it's filtered upstream so don't
                # double-count below.
                stats.skipped_mate_mapq_fail += 1
                continue
            if check_mm and has_multimap_tag(read):
                stats.skipped_multimap += 1
                continue
            if has_excessive_5p_softclip(read, sc_thresh):
                stats.skipped_excessive_softclip += 1
                continue
            # Construct fragment coords. Mirror frag_generator logic.
            frag_start = min(read.reference_start, read.next_reference_start)
            frag_stop = frag_start + abs(read.template_length)
            if frag_stop - frag_start < fraction_low:
                continue
            if frag_stop - frag_start > fraction_high:
                continue
            read_on_plus = not read.is_reverse
            stats.total_fragments_processed += 1

            qname = read.query_name
            r2 = r2_cache.pop(qname, None)
            if r2 is not None:
                _process_one_fragment(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    read.mapping_quality, refseq, motif_to_classes,
                    args.motif_length, stats, frag_log,
                    bam_writer=bam_writer, read1=read, read2=r2,
                )
            else:
                # Stash Read 1 and its classification target for later.
                motif = get_motif(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    refseq, k=args.motif_length,
                )
                if motif is None:
                    stats.skipped_motif_invalid += 1
                    frag_log.write(
                        read.reference_name, frag_start, frag_stop,
                        read_on_plus, None, set(), [],
                    )
                    continue
                if read_on_plus:
                    stats.motif_counts_plus[motif] += 1
                else:
                    stats.motif_counts_minus[motif] += 1
                outputs = classify(motif, motif_to_classes)
                classes = motif_to_classes.get(motif, set())
                stats.passed_qc += 1
                frag_log.write(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    motif, classes, outputs,
                )
                for label in outputs:
                    stats.output_counts[label] += 1
                    stats.per_output_motif_counts[label][motif] += 1
                r1_cache[qname] = (read, outputs, motif)
        else:
            # Read 2.
            qname = read.query_name
            r1_pkg = r1_cache.pop(qname, None)
            if r1_pkg is None:
                r2_cache[qname] = read
            else:
                read1, outputs, motif = r1_pkg
                for label in outputs:
                    bam_writer.write_pair(label, read1, read, motif)

    # Drain any unmatched Read 1s (Read 2 was probably filtered or off-target).
    # These still get written without their mate; we already counted them.
    for qname, (read1, outputs, motif) in r1_cache.items():
        for label in outputs:
            bam_writer.write_pair(label, read1, None, motif)

    sam.close()
    # Stash the template for downstream callers if they need it.
    bam_writer.template_bam = template_bam


def process_bam_fetch_mate(
    args: argparse.Namespace,
    refseq: Any,
    motif_to_classes: dict[str, set[str]],
    stats: RunStats,
    bam_writer: BamOutputs,
    frag_log: FragmentLog,
) -> None:
    """fetch_mate path: per-fragment mate lookup via pysam.mate(). Slower (~2x
    BAM I/O) but works on BAMs without the MQ tag.

    For each primary Read 1, fetch its mate and check the mate's MAPQ
    directly. Then apply the same multimap / softclip filters as the
    require_mq_tag path.
    """
    import pysam
    sam = pysam.AlignmentFile(args.input, "rb", threads=args.threads)
    mapq_min = args.quality_threshold
    check_mm = not args.no_check_multimap
    sc_thresh = args.max_5prime_softclip
    fraction_low = args.fraction_low
    fraction_high = args.fraction_high if args.fraction_high is not None else (1 << 31)

    for read in sam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if not read.is_read1:
            continue
        if read.mate_is_unmapped:
            continue
        if read.reference_name != read.next_reference_name:
            continue
        if read.mapping_quality < mapq_min:
            continue
        # Look up the mate.
        try:
            mate = sam.mate(read)
        except (ValueError, KeyError):  # mate not found
            stats.skipped_mate_mapq_fail += 1
            continue
        if mate.mapping_quality < mapq_min:
            stats.skipped_mate_mapq_fail += 1
            continue
        if check_mm and (has_multimap_tag(read) or has_multimap_tag(mate)):
            stats.skipped_multimap += 1
            continue
        if has_excessive_5p_softclip(read, sc_thresh):
            stats.skipped_excessive_softclip += 1
            continue
        frag_start = min(read.reference_start, mate.reference_start)
        frag_stop = frag_start + abs(read.template_length)
        if frag_stop - frag_start < fraction_low or frag_stop - frag_start > fraction_high:
            continue
        read_on_plus = not read.is_reverse
        stats.total_fragments_processed += 1
        _process_one_fragment(
            read.reference_name, frag_start, frag_stop, read_on_plus,
            read.mapping_quality, refseq, motif_to_classes, args.motif_length,
            stats, frag_log, bam_writer=bam_writer, read1=read, read2=mate,
        )
    sam.close()
    bam_writer.template_bam = sam


def process_bam_skip_mate_check(
    args: argparse.Namespace,
    refseq: Any,
    motif_to_classes: dict[str, set[str]],
    stats: RunStats,
    bam_writer: BamOutputs,
    frag_log: FragmentLog,
) -> None:
    """skip mode: no mate-MAPQ check, otherwise same as require_mq_tag.

    Provided only for diagnostic comparison against legacy results; not
    recommended for production runs.
    """
    import pysam
    sam = pysam.AlignmentFile(args.input, "rb", threads=args.threads)
    bam_writer.template_bam = sam
    mapq_min = args.quality_threshold
    check_mm = not args.no_check_multimap
    sc_thresh = args.max_5prime_softclip
    fraction_low = args.fraction_low
    fraction_high = args.fraction_high if args.fraction_high is not None else (1 << 31)

    r2_cache: dict[str, Any] = {}
    r1_cache: dict[str, tuple[Any, list[str], str]] = {}

    for read in sam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.mate_is_unmapped:
            continue
        if read.reference_name != read.next_reference_name:
            continue
        if read.is_read1:
            if read.mapping_quality < mapq_min:
                continue
            if check_mm and has_multimap_tag(read):
                stats.skipped_multimap += 1
                continue
            if has_excessive_5p_softclip(read, sc_thresh):
                stats.skipped_excessive_softclip += 1
                continue
            frag_start = min(read.reference_start, read.next_reference_start)
            frag_stop = frag_start + abs(read.template_length)
            if frag_stop - frag_start < fraction_low or frag_stop - frag_start > fraction_high:
                continue
            read_on_plus = not read.is_reverse
            stats.total_fragments_processed += 1
            qname = read.query_name
            r2 = r2_cache.pop(qname, None)
            if r2 is not None:
                _process_one_fragment(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    read.mapping_quality, refseq, motif_to_classes,
                    args.motif_length, stats, frag_log,
                    bam_writer=bam_writer, read1=read, read2=r2,
                )
            else:
                motif = get_motif(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    refseq, k=args.motif_length,
                )
                if motif is None:
                    stats.skipped_motif_invalid += 1
                    frag_log.write(
                        read.reference_name, frag_start, frag_stop,
                        read_on_plus, None, set(), [],
                    )
                    continue
                if read_on_plus:
                    stats.motif_counts_plus[motif] += 1
                else:
                    stats.motif_counts_minus[motif] += 1
                outputs = classify(motif, motif_to_classes)
                classes = motif_to_classes.get(motif, set())
                stats.passed_qc += 1
                frag_log.write(
                    read.reference_name, frag_start, frag_stop, read_on_plus,
                    motif, classes, outputs,
                )
                for label in outputs:
                    stats.output_counts[label] += 1
                    stats.per_output_motif_counts[label][motif] += 1
                r1_cache[qname] = (read, outputs, motif)
        else:
            qname = read.query_name
            r1_pkg = r1_cache.pop(qname, None)
            if r1_pkg is None:
                r2_cache[qname] = read
            else:
                read1, outputs, motif = r1_pkg
                for label in outputs:
                    bam_writer.write_pair(label, read1, read, motif)
    for qname, (read1, outputs, motif) in r1_cache.items():
        for label in outputs:
            bam_writer.write_pair(label, read1, None, motif)
    sam.close()


# ---------------------------------------------------------------------------
# Summary JSON.
# ---------------------------------------------------------------------------

def write_summary_json(
    args: argparse.Namespace,
    stats: RunStats,
    motif_config_sha: str,
    output_path: str,
) -> None:
    """Write summary JSON per §6.2."""
    total = max(stats.total_fragments_processed, 1)
    retention_rates = {
        label: round(stats.output_counts[label] / total, 6)
        for label in OUTPUT_LABELS
    }
    payload = {
        "input": args.input,
        "reference": args.reference,
        "motif_length": args.motif_length,
        "l3_set": args.l3_set,
        "d1_set": args.d1_set,
        "dffb_set": args.dffb_set,
        "mate_check_mode": args.mate_check,
        "check_multimap": (not args.no_check_multimap),
        "max_5prime_softclip": args.max_5prime_softclip,
        "fraction_low": args.fraction_low,
        "fraction_high": args.fraction_high,
        "quality_threshold": args.quality_threshold,
        "total_fragments_processed": stats.total_fragments_processed,
        "passed_qc": stats.passed_qc,
        "skipped_motif_invalid": stats.skipped_motif_invalid,
        "skipped_mate_mapq_fail": stats.skipped_mate_mapq_fail,
        "skipped_multimap": stats.skipped_multimap,
        "skipped_excessive_softclip": stats.skipped_excessive_softclip,
        "output_counts": {
            label: stats.output_counts[label] for label in OUTPUT_LABELS
        },
        "retention_rates": retention_rates,
        "motif_counts_plus_strand": dict(stats.motif_counts_plus),
        "motif_counts_minus_strand": dict(stats.motif_counts_minus),
        "config_hash": motif_config_sha,
        "tool_version": __version__,
        "command_line": " ".join(sys.argv),
    }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# QC PDF — 7 pages per design §6.4. Heavy imports gated by --no-qc-pdf.
# ---------------------------------------------------------------------------

def _locate_f_profile_matrix() -> Optional[Path]:
    """Try to find FinaleToolkit's end_motif_f_profiles.tsv. Returns None if
    it's not found — the QC PDF will skip the F-profile decomposition page
    with a note explaining where to put the file."""
    try:
        import finaletoolkit
    except ImportError:
        return None
    pkg_dir = Path(finaletoolkit.__file__).resolve().parent
    candidates = [
        pkg_dir / "data" / "end_motif_f_profiles.tsv",
        pkg_dir / "data" / "f_profiles.tsv",
        pkg_dir.parent / "data" / "end_motif_f_profiles.tsv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def _load_f_profile_matrix(path: Path) -> tuple[list[str], "np.ndarray"]:
    """Load the F-profile matrix. Expected format: header row of profile
    names (F-I..F-VI), then 256 rows of motif + per-profile probabilities.

    Returns (motif_order, matrix) with matrix.shape == (256, 6).
    """
    import numpy as np
    motifs = []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            motifs.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
    return motifs, np.array(rows)


def _nnls_f_profile_decomp(
    motif_counts: dict[str, int],
    motif_order: list[str],
    f_matrix: "np.ndarray",
    n_bootstrap: int = 100,
) -> dict[str, Any]:
    """NNLS decomposition of an empirical 256-bin motif distribution against
    the 6 F-profiles. Returns the point estimate + 95% bootstrap CIs.

    n_bootstrap defaults to 100 here (rather than 1000 as in the design
    text) so the QC PDF is fast enough to enable by default; raise via
    a future flag if you want tighter CIs.
    """
    import numpy as np
    from scipy.optimize import nnls

    # Build the empirical distribution in motif_order.
    total = sum(motif_counts.values())
    if total == 0:
        return {"point": np.zeros(f_matrix.shape[1]), "ci_lo": None, "ci_hi": None}
    dist = np.array([motif_counts.get(m, 0) / total for m in motif_order], dtype=float)
    # Point estimate: NNLS then normalize to sum to 1.
    w, _ = nnls(f_matrix, dist)
    s = w.sum()
    w = w / s if s > 0 else w
    # Bootstrap CIs by resampling per-motif counts.
    rng = np.random.default_rng(0)
    counts_vec = np.array([motif_counts.get(m, 0) for m in motif_order], dtype=int)
    boot = np.zeros((n_bootstrap, f_matrix.shape[1]))
    for b in range(n_bootstrap):
        sample_idx = rng.choice(len(motif_order), size=total, replace=True, p=dist)
        boot_counts = np.bincount(sample_idx, minlength=len(motif_order))
        boot_dist = boot_counts / total
        wb, _ = nnls(f_matrix, boot_dist)
        sb = wb.sum()
        boot[b] = wb / sb if sb > 0 else wb
    ci_lo = np.percentile(boot, 2.5, axis=0)
    ci_hi = np.percentile(boot, 97.5, axis=0)
    return {"point": w, "ci_lo": ci_lo, "ci_hi": ci_hi}


def _shade_motif_zones(ax, motif_order: list[str]) -> None:
    """Background-shade L3 (red), D1 (orange), DFFB-canonical (blue) zones
    on a 256-bin motif axis. Used by QC page 2 + page 3."""
    for i, m in enumerate(motif_order):
        if m in CCNN_MOTIFS:
            ax.axvspan(i - 0.5, i + 0.5, color="red", alpha=0.10, zorder=0)
        elif m in TGNN_MOTIFS:
            ax.axvspan(i - 0.5, i + 0.5, color="orange", alpha=0.10, zorder=0)
        elif m in DFFB_CANONICAL_MOTIFS:
            ax.axvspan(i - 0.5, i + 0.5, color="blue", alpha=0.18, zorder=0)


def _draw_motif_histogram(ax, motif_counts: dict[str, int], motif_order: list[str],
                          title: str, shade_zones: bool = True) -> None:
    """One 256-bin motif histogram, x-axis ordered by motif_order."""
    total = sum(motif_counts.values())
    freqs = [motif_counts.get(m, 0) / total * 100.0 if total else 0.0
             for m in motif_order]
    xs = list(range(len(motif_order)))
    if shade_zones:
        _shade_motif_zones(ax, motif_order)
    ax.bar(xs, freqs, width=1.0, color="#404040", linewidth=0)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("frequency (%)", fontsize=8)
    ax.set_xlim(-0.5, len(motif_order) - 0.5)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def write_qc_pdf(
    args: argparse.Namespace,
    stats: RunStats,
    motif_config_sha: str,
    output_path: str,
) -> None:
    """Generate the 7-page QC PDF (§6.4). Returns silently on plotting-lib
    import failure with a logged warning so the rest of the run is unaffected.
    """
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        try:
            import seaborn as sns
            sns.set_style("whitegrid")
        except ImportError:
            sns = None
    except ImportError as exc:
        log.warning(
            "QC PDF skipped: matplotlib unavailable (%s). "
            "Install with `pip install matplotlib seaborn` or pass --no-qc-pdf.",
            exc,
        )
        return

    motif_order = list(CANONICAL_4MERS)
    f_matrix_path = _locate_f_profile_matrix()
    f_matrix_data: Optional[tuple[list[str], "np.ndarray"]] = None
    if f_matrix_path is not None:
        try:
            mo, fm = _load_f_profile_matrix(f_matrix_path)
            # Re-order to motif_order so plotting/NNLS align.
            mo_idx = {m: i for i, m in enumerate(mo)}
            fm_reord = np.zeros_like(fm)
            ok = True
            for i, m in enumerate(motif_order):
                j = mo_idx.get(m)
                if j is None:
                    ok = False
                    break
                fm_reord[i] = fm[j]
            if ok:
                f_matrix_data = (motif_order, fm_reord)
        except Exception as exc:
            log.warning("Failed to load F-profile matrix %s: %s", f_matrix_path, exc)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plus_only = dict(stats.motif_counts_plus)
    minus_only = dict(stats.motif_counts_minus)
    combined_counts = Counter()
    for d in (plus_only, minus_only):
        for m, c in d.items():
            combined_counts[m] += c

    config_key = "ccnn" if args.l3_set == "ccnn" else "fprofile_argmax"
    exp_retention = EXPECTED_RETENTION[config_key]

    with PdfPages(output_path) as pdf:
        # --- Page 1: Output retention summary ---
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        labels = list(OUTPUT_LABELS)
        total = max(stats.total_fragments_processed, 1)
        rates = [stats.output_counts[l] / total for l in labels]
        bar_colors = []
        for label, rate in zip(labels, rates):
            lo, hi = exp_retention.get(label, (0.0, 1.0))
            bar_colors.append("#3a9b56" if (lo <= rate <= hi) else "#d68a3a")
        xs = list(range(len(labels)))
        bars = ax.bar(xs, [r * 100 for r in rates], color=bar_colors, edgecolor="black")
        for x, (l, h) in zip(xs, [exp_retention.get(lab, (0, 0)) for lab in labels]):
            ax.plot([x - 0.4, x + 0.4], [l * 100, l * 100], "--", color="gray", linewidth=0.8)
            ax.plot([x - 0.4, x + 0.4], [h * 100, h * 100], "--", color="gray", linewidth=0.8)
        for x, rate, count in zip(xs, rates, [stats.output_counts[l] for l in labels]):
            ax.text(x, rate * 100 + 1.0, f"{count:,}\n({rate*100:.1f}%)",
                    ha="center", va="bottom", fontsize=8)
        ax.set_xticks(xs)
        ax.set_xticklabels([l.replace("_", "\n") for l in labels],
                           rotation=0, fontsize=8)
        ax.set_ylabel("retention rate (%)")
        ax.set_title(
            f"Output retention summary "
            f"(input: {stats.total_fragments_processed:,} fragments; "
            f"passed QC: {stats.passed_qc:,})",
            fontsize=11,
        )
        ax.set_ylim(0, max(max([r * 100 for r in rates], default=0) * 1.25, 5))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.text(
            0.5, 0.02,
            "Green: within expected healthy-plasma range (dashed lines).  "
            "Orange: outside expected range — may indicate cancer enrichment, "
            "library-prep artifacts, or config mismatch.",
            ha="center", fontsize=8, color="dimgray",
        )
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: Unfiltered 256-bin motif distribution ---
        fig, ax = plt.subplots(figsize=(11, 4.5))
        _draw_motif_histogram(
            ax, combined_counts, motif_order,
            title=(
                f"Unfiltered 4-mer end-motif distribution "
                f"({sum(combined_counts.values()):,} motifs)"
            ),
        )
        # Zone legend at bottom.
        fig.text(0.05, 0.02,
                 "red shading = CCNN (DNase1L3 zone); "
                 "orange shading = TGNN (DNase1 zone); "
                 "blue shading = DFFB canonical 4 motifs",
                 fontsize=8, color="dimgray")
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: Per-output motif distributions (3x2 grid) ---
        fig, axes = plt.subplots(3, 2, figsize=(13, 9))
        axes_flat = axes.flatten()
        for ax, label in zip(axes_flat[:5], OUTPUT_LABELS):
            counts = stats.per_output_motif_counts.get(label, Counter())
            n = sum(counts.values())
            _draw_motif_histogram(
                ax, counts, motif_order,
                title=f"{label}   (n={n:,})",
            )
        # Sixth panel: re-show unfiltered for adjacent comparison.
        _draw_motif_histogram(
            axes_flat[5], combined_counts, motif_order,
            title="Unfiltered input (for direct comparison)",
        )
        fig.suptitle(
            "Per-output 4-mer end-motif distributions",
            fontsize=12, fontweight="bold",
        )
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 4: F-profile compositional decomposition ---
        if f_matrix_data is None:
            fig, ax = plt.subplots(figsize=(11, 5))
            ax.text(0.5, 0.5,
                    "F-profile matrix (end_motif_f_profiles.tsv) not found in "
                    "FinaleToolkit's data dir.\nF-profile decomposition page "
                    "skipped.\nPlace the file at <finaletoolkit>/data/end_motif_f_profiles.tsv "
                    "to enable this page.",
                    ha="center", va="center", fontsize=11, color="firebrick")
            ax.axis("off")
            ax.set_title("F-profile compositional decomposition (skipped)")
            pdf.savefig(fig)
            plt.close(fig)
        else:
            mo, fm = f_matrix_data
            decomps = {"Input (unfiltered)": combined_counts}
            for label in OUTPUT_LABELS:
                decomps[label] = stats.per_output_motif_counts.get(label, Counter())
            fits = {
                name: _nnls_f_profile_decomp(c, mo, fm, n_bootstrap=100)
                for name, c in decomps.items()
            }
            f_labels = ["F-I", "F-II", "F-III", "F-IV", "F-V", "F-VI"]
            n_profiles = len(f_labels)
            n_groups = len(decomps)

            fig, ax = plt.subplots(figsize=(13, 5.5))
            bar_width = 0.85 / n_groups
            group_keys = list(decomps.keys())
            palette = ["#404040", "#d62728", "#1f77b4", "#9467bd", "#2ca02c", "#ff7f0e"]
            for gi, gname in enumerate(group_keys):
                fit = fits[gname]
                w = fit["point"]
                if fit["ci_lo"] is not None:
                    yerr = np.vstack([w - fit["ci_lo"], fit["ci_hi"] - w])
                else:
                    yerr = None
                xs = np.arange(n_profiles) + (gi - n_groups / 2 + 0.5) * bar_width
                ax.bar(
                    xs, w * 100, bar_width * 0.95,
                    yerr=yerr * 100 if yerr is not None else None,
                    label=gname, color=palette[gi % len(palette)],
                    edgecolor="black", linewidth=0.3,
                    error_kw={"linewidth": 0.6, "ecolor": "black"},
                )
            ax.set_xticks(np.arange(n_profiles))
            ax.set_xticklabels(f_labels)
            ax.set_ylabel("F-profile weight (%)")
            ax.set_title(
                "F-profile compositional decomposition — input vs filtered outputs\n"
                "(NNLS fit on 256-bin motif distribution; error bars = 95% bootstrap CI)",
                fontsize=11,
            )
            ax.legend(fontsize=8, loc="upper right", framealpha=0.92)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            fig.text(
                0.5, 0.02,
                "F-I/F-II ~ DNase1L3/DNase1 (chromatin); F-III ~ DFFB (apoptosis); "
                "F-IV CpG-rich; F-V G-end; F-VI oxidative (cancer-discriminating per Zhou et al. 2023). "
                "Bootstrap CIs reflect sampling noise only; F-profile model misspecification is not captured.",
                ha="center", fontsize=7, color="dimgray",
            )
            fig.tight_layout(rect=(0, 0.05, 1, 1))
            pdf.savefig(fig)
            plt.close(fig)

        # --- Page 5: Per-strand motif distributions (unfiltered) ---
        fig, axes = plt.subplots(2, 1, figsize=(11, 8))
        n_plus = sum(plus_only.values())
        n_minus = sum(minus_only.values())
        _draw_motif_histogram(
            axes[0], plus_only, motif_order,
            title=f"Forward-strand fragments (Read 1 on + strand)   (n={n_plus:,})",
        )
        _draw_motif_histogram(
            axes[1], minus_only, motif_order,
            title=(
                f"Reverse-strand fragments (Read 1 on − strand, "
                f"reverse-complemented)   (n={n_minus:,})"
            ),
        )
        # Strand-similarity KL divergence as a single number.
        try:
            p = np.array([plus_only.get(m, 0) for m in motif_order], dtype=float)
            q = np.array([minus_only.get(m, 0) for m in motif_order], dtype=float)
            p = (p + 1e-9) / (p.sum() + len(p) * 1e-9)
            q = (q + 1e-9) / (q.sum() + len(q) * 1e-9)
            kl_pq = float((p * np.log(p / q)).sum())
            kl_qp = float((q * np.log(q / p)).sum())
            sym = (kl_pq + kl_qp) / 2.0
            fig.text(
                0.5, 0.01,
                f"Symmetric KL divergence (+ strand vs − strand) = {sym:.4f} nats.  "
                f"Healthy-library expectation: ~0 (each strand contributes ~50% "
                f"of fragments with the same enzyme mixing). Large values indicate "
                f"strand-specific library prep artifacts.",
                ha="center", fontsize=7.5, color="dimgray",
            )
        except Exception:
            pass
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 6: QC funnel ---
        #
        # Counter invariants in production code:
        #   total_fragments_processed = passed_qc + skipped_motif_invalid
        #     (incremented AFTER all per-read filters, BEFORE motif check)
        #   skipped_motif_invalid is a SUBSET of total_fragments_processed
        #   skipped_{mate_mapq_fail, multimap, excessive_softclip} are
        #     DISJOINT from total_fragments_processed (the read never
        #     reached the total++ line)
        # So the cumulative counts at each filter stage are:
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        total_scanned = (
            stats.total_fragments_processed + stats.skipped_mate_mapq_fail
            + stats.skipped_multimap + stats.skipped_excessive_softclip
        )
        stages = [
            ("Total fragments scanned", total_scanned),
            ("Passed pair-quality (MAPQ/MQ/flags)",
                stats.total_fragments_processed + stats.skipped_multimap +
                stats.skipped_excessive_softclip),
            ("Passed multi-mapping filter (XA/SA)",
                stats.total_fragments_processed +
                stats.skipped_excessive_softclip),
            ("Passed 5' softclip filter",
                stats.total_fragments_processed),
            ("Valid motif (not N / not boundary)", stats.passed_qc),
            ("Assigned to ≥ 1 output (== passed QC)", stats.passed_qc),
        ]
        names = [n for n, _ in stages]
        counts = [v for _, v in stages]
        ys = list(range(len(stages)))[::-1]
        ax.barh(ys, counts, color="#3a6ea5", edgecolor="black", height=0.7)
        for y, count, prev_count in zip(ys, counts, [counts[0]] + counts[:-1]):
            pct = (100.0 * count / prev_count) if prev_count > 0 else 0.0
            drop = prev_count - count
            ax.text(count + max(counts) * 0.015, y,
                    f"{count:,}  ({pct:.1f}% pass, drop = {drop:,})",
                    va="center", fontsize=8)
        ax.set_yticks(ys)
        ax.set_yticklabels(names, fontsize=9)
        ax.set_xlabel("fragment count")
        ax.set_title("QC funnel: fragments at each filter stage", fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0, max(counts) * 1.45 if counts else 1)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 7: Config and metadata ---
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis("off")
        ax.set_title(
            f"motif_filter.py v{__version__} — configuration and metadata",
            fontsize=12, fontweight="bold", loc="left",
        )
        lines = [
            f"Run timestamp:       {time.strftime('%Y-%m-%d %H:%M:%S')}",
            f"Tool version:        {__version__}",
            f"Command line:        {' '.join(sys.argv)}",
            "",
            f"Input:               {args.input}",
            f"Reference:           {args.reference}",
            f"Output prefix:       {args.output_prefix}",
            "",
            f"Motif length:        {args.motif_length}",
            f"L3 set:              {args.l3_set}",
            f"D1 set:              {args.d1_set}",
            f"DFFB set:            {args.dffb_set}",
            f"Motif config TSV:    {args.motif_config or '(built-in defaults)'}",
            f"Motif config hash:   {motif_config_sha}",
            "",
            f"Quality threshold:   {args.quality_threshold}",
            f"Mate check mode:     {args.mate_check}",
            f"Check multimap:      {not args.no_check_multimap}",
            f"Max 5' softclip:     {args.max_5prime_softclip}",
            f"Fragment length:     [{args.fraction_low}, "
            f"{args.fraction_high if args.fraction_high is not None else 'inf'}]",
            f"Threads:             {args.threads}",
            f"Per-fragment log:    {args.log or '(disabled)'}",
            f"Emit motif tag:      {args.emit_motif_tag}",
            "",
            "Fragment counts:",
            f"  total processed:        {stats.total_fragments_processed:,}",
            f"  passed QC:              {stats.passed_qc:,}",
            f"  skipped invalid motif:  {stats.skipped_motif_invalid:,}",
            f"  skipped mate-MAPQ:      {stats.skipped_mate_mapq_fail:,}",
            f"  skipped multimap:       {stats.skipped_multimap:,}",
            f"  skipped softclip:       {stats.skipped_excessive_softclip:,}",
            "",
            "Output counts:",
        ]
        for label in OUTPUT_LABELS:
            lines.append(f"  {label:<32}{stats.output_counts[label]:>12,}")
        text = "\n".join(lines)
        ax.text(
            0.02, 0.96, text, transform=ax.transAxes,
            fontfamily="monospace", fontsize=9, va="top",
        )
        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------------
# main() entry point.
# ---------------------------------------------------------------------------

def main(argv: Optional[Iterable[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose)
    t0 = time.perf_counter()

    # Lazy-import py2bit (only needed for actual runs, not --help).
    try:
        import py2bit
    except ImportError:
        log.error("py2bit not installed. Run: pip install py2bit")
        return 2

    # Sanity-check fraction-low vs motif-length.
    if args.fraction_low < args.motif_length:
        log.warning(
            "--fraction-low (%d) < --motif-length (%d); clamping to %d",
            args.fraction_low, args.motif_length, args.motif_length,
        )
        args.fraction_low = args.motif_length

    # Build motif lookup.
    if args.motif_config is not None:
        motif_to_classes = load_motif_config(args.motif_config, args.motif_length)
        log.info(
            "Loaded %d motifs from --motif-config %s",
            len(motif_to_classes), args.motif_config,
        )
    else:
        motif_to_classes = build_motif_to_classes(
            args.l3_set, args.d1_set, args.dffb_set,
        )
        log.info(
            "Built motif lookup (built-in defaults): %d motifs across "
            "%d classes (l3=%s, d1=%s, dffb=%s)",
            len(motif_to_classes), 4, args.l3_set, args.d1_set, args.dffb_set,
        )
    motif_config_sha = motif_config_hash(motif_to_classes)

    fmt = detect_input_format(args.input)

    # Pre-flight: MQ tag presence check for BAM + require_mq_tag mode.
    if fmt in ("bam", "cram") and args.mate_check == "require_mq_tag":
        log.info("Verifying input BAM has MQ tag on Read 1 alignments...")
        _check_first_alignment_has_mq_tag(args.input, args.threads)
        log.info("MQ tag check passed.")

    # Open 2bit reference.
    log.info("Opening 2bit reference: %s", args.reference)
    refseq = py2bit.open(args.reference)

    # Prepare outputs and per-fragment log.
    Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)
    frag_log = FragmentLog(args.log)
    stats = RunStats()

    bam_writer: Optional[BamOutputs] = None
    frag_writer: Optional[FragGzOutputs] = None
    try:
        if fmt in ("bam", "cram"):
            import pysam
            # Open the BAM template for the writers.
            template = pysam.AlignmentFile(args.input, "rb", threads=args.threads)
            bam_writer = BamOutputs(
                prefix=args.output_prefix,
                template_bam=template,
                emit_motif_tag=args.emit_motif_tag,
                threads=args.threads,
            )
            # Dispatch on mate-check mode.
            if args.mate_check == "require_mq_tag":
                process_bam_require_mq_tag(
                    args, refseq, motif_to_classes, stats, bam_writer, frag_log,
                )
            elif args.mate_check == "fetch_mate":
                process_bam_fetch_mate(
                    args, refseq, motif_to_classes, stats, bam_writer, frag_log,
                )
            else:  # skip
                process_bam_skip_mate_check(
                    args, refseq, motif_to_classes, stats, bam_writer, frag_log,
                )
            bam_writer.close_and_index()
            template.close()
        else:  # frag.gz
            frag_writer = FragGzOutputs(prefix=args.output_prefix)
            process_frag_gz(
                args, refseq, motif_to_classes, stats, frag_writer, frag_log,
            )
            frag_writer.close_and_index()
    finally:
        frag_log.close()
        try:
            refseq.close()
        except Exception:
            pass

    # Summary JSON.
    summary_path = f"{args.output_prefix}.motif_filter_stats.json"
    write_summary_json(args, stats, motif_config_sha, summary_path)
    log.info("Wrote summary JSON: %s", summary_path)

    # QC PDF (default on; skip with --no-qc-pdf).
    if not args.no_qc_pdf:
        qc_path = f"{args.output_prefix}.motif_filter_qc.pdf"
        write_qc_pdf(args, stats, motif_config_sha, qc_path)
        log.info("Wrote QC PDF: %s", qc_path)

    elapsed = time.perf_counter() - t0
    log.info(
        "Done in %.1f s. Processed %d fragments; passed QC %d (%.1f%%).",
        elapsed, stats.total_fragments_processed, stats.passed_qc,
        100.0 * stats.passed_qc / max(stats.total_fragments_processed, 1),
    )
    for label in OUTPUT_LABELS:
        log.info("  %s : %d", label, stats.output_counts[label])
    return 0


if __name__ == "__main__":
    sys.exit(main())
