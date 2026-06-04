#!/usr/bin/env python3
"""Build SIVmac239m2-based synthetic test fixtures for umi-pipeline-nf-HIV."""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys

import pysam

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIVMAC_DIR = os.path.join(REPO_ROOT, "tests", "input", "sivmac239m2")
REF_FA = os.path.join(SIVMAC_DIR, "sivmac239m2.ltr_masked.fa")
CONTIG = "SIVmac239m2"
LN_START, LN_END = 1179, 1774
SLN_START, SLN_END = 251, 6900
DELETION_BP = 250

FWD_UMI = "ACGTACGTACGTACGTAC"
REV_UMI = "TGCATGCATGCATGCATG"
FWD_CONTEXT = "CAAGCAGAAGACGGCATACGAGAT"
REV_CONTEXT = "TAGGGAGCCGTCAGGATCAG"


def fetch_region(ref_path: str, start: int, end: int) -> str:
    with pysam.FastaFile(ref_path) as ff:
        return ff.fetch(CONTIG, start, end)


def apply_deletion(seq: str, del_bp: int, del_start_frac: float = 0.35) -> str:
    del_start = int(len(seq) * del_start_frac)
    del_end = min(del_start + del_bp, len(seq))
    return seq[:del_start] + seq[del_end:]


def combine_umi_tag(strand: str) -> str:
    if strand == "+":
        return FWD_UMI + REV_UMI
    return _rev_comp(REV_UMI) + _rev_comp(FWD_UMI)


def _rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def write_cluster_fasta(path: str, reads: list[dict]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        for i, r in enumerate(reads):
            name = (
                f"read{i};strand={r['strand']};umi_fwd_dist=0;umi_rev_dist=0;"
                f"umi_fwd_seq={FWD_UMI};umi_rev_seq={REV_UMI};seq={r['seq']}"
            )
            fh.write(f">{name}\n{r['tag']}\n")


def build_length_cluster_fixtures(out_dir: str, ln_seq: str, ln_del: str) -> None:
    tag = combine_umi_tag("+")
    full_reads = [{"strand": "+", "seq": ln_seq, "tag": tag} for _ in range(6)]
    del_reads = [{"strand": "+", "seq": ln_del, "tag": tag} for _ in range(6)]
    near_reads = [
        {"strand": "+", "seq": ln_seq[: len(ln_seq) - 20] + "A" * 20, "tag": tag}
        for _ in range(4)
    ]
    write_cluster_fasta(
        os.path.join(out_dir, "cluster_deletion_pair"),
        full_reads + del_reads,
    )
    write_cluster_fasta(
        os.path.join(out_dir, "cluster_within_threshold"),
        full_reads[:4] + near_reads,
    )


def build_split_reads_bam(out_dir: str, ref_path: str, ln_seq: str, ln_del: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    bed_path = os.path.join(out_dir, "amplicon.bed")
    with open(bed_path, "w") as fh:
        fh.write(f"{CONTIG}\t{LN_START}\t{LN_END}\tamplicon\n")

    bam_path = os.path.join(out_dir, "merged.bam")
    if os.path.exists(bam_path):
        os.remove(bam_path)

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": len(fetch_region(ref_path, 0, 10_000)), "SN": CONTIG}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=pysam.AlignmentHeader.from_dict(header)) as bam:
        _add_deletion_mapped_read(
            bam, name="del_0", seq=ln_del, ref_start=LN_START, del_ref_len=DELETION_BP
        )
        _add_mapped_read(bam, name="full_0", seq=ln_seq, ref_start=LN_START, is_reverse=False)
        for i in range(1, 4):
            _add_mapped_read(
                bam, name=f"full_{i}", seq=ln_seq, ref_start=LN_START + i * 3, is_reverse=False
            )

    pysam.index(bam_path)
    shutil.copy(bed_path, os.path.join(SIVMAC_DIR, "sivmac239m2_ln.bed"))


def _add_mapped_read(bam, name: str, seq: str, ref_start: int, is_reverse: bool) -> None:
    a = pysam.AlignedSegment(bam.header)
    a.query_name = name
    a.query_sequence = seq
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigar = [(0, len(seq))]
    a.is_reverse = is_reverse
    a.is_secondary = False
    a.is_supplementary = False
    bam.write(a)


def _add_deletion_mapped_read(
    bam, name: str, seq: str, ref_start: int, del_ref_len: int
) -> None:
    """Query with internal deletion; reference alignment still spans the BED."""
    match_left = 80
    match_right = len(seq) - match_left
    a = pysam.AlignedSegment(bam.header)
    a.query_name = name
    a.query_sequence = seq
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigar = [(0, match_left), (2, del_ref_len), (0, match_right)]
    a.is_reverse = False
    a.is_secondary = False
    a.is_supplementary = False
    bam.write(a)


def build_smoke_fastq(out_dir: str, ln_seq: str, ln_del: str) -> None:
    """Synthetic nanopore reads with UMI flanks for short LN amplicon smoke test."""
    fastq_dir = os.path.join(out_dir, "fastq_pass", "barcode01")
    os.makedirs(fastq_dir, exist_ok=True)
    fastq_path = os.path.join(fastq_dir, "sivmac_smoke.fastq")

    rng = random.Random(42)
    records = []
    for i in range(25):
        inner = ln_seq if i % 2 == 0 else ln_del
        seq = FWD_CONTEXT + FWD_UMI + inner + REV_UMI + REV_CONTEXT
        qual = "".join(chr(rng.randint(33, 70)) for _ in seq)
        records.append((f"smoke_{i}", seq, qual))

    with open(fastq_path, "w") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", default=REF_FA)
    args = parser.parse_args()

    if not os.path.isfile(args.ref):
        print(f"Missing reference: {args.ref}", file=sys.stderr)
        return 1

    ln_seq = fetch_region(args.ref, LN_START, LN_END)
    ln_del = apply_deletion(ln_seq, DELETION_BP)

    length_dir = os.path.join(
        REPO_ROOT, "tests", "input", "reformat_filter_cluster", "sivmac", "length_cluster"
    )
    build_length_cluster_fixtures(length_dir, ln_seq, ln_del)

    split_dir = os.path.join(REPO_ROOT, "tests", "input", "split_reads", "sivmac_deletion")
    build_split_reads_bam(split_dir, args.ref, ln_seq, ln_del)

    smoke_dir = os.path.join(REPO_ROOT, "tests", "input", "pipeline", "sivmac_smoke")
    build_smoke_fastq(smoke_dir, ln_seq, ln_del)

    print("Wrote fixtures:")
    print(f"  {length_dir}")
    print(f"  {split_dir}")
    print(f"  {smoke_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
