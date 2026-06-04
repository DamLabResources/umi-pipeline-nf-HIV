#!/usr/bin/env python3
"""Unit tests for SIVmac split-read fixtures (no Nextflow)."""

import glob
import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
from filter_reads import parse_args, filter_reads, is_long_read, is_short_read  # noqa: E402

BED = os.path.join(REPO, "tests/input/split_reads/sivmac_deletion/amplicon.bed")
BAM = os.path.join(REPO, "tests/input/split_reads/sivmac_deletion/merged.bam")


class TestSivmacSplitReads(unittest.TestCase):
    def _run(self, mode):
        tmp = tempfile.mkdtemp()
        args = parse_args(
            [
                "-o",
                tmp,
                "--min_overlap",
                "0.95",
                "--adapter_length",
                "100",
                "--split_read_filter_mode",
                mode,
                "--output_format",
                "fasta",
                "--output_filename",
                "t",
                BED,
                BAM,
            ]
        )
        filter_reads(args)
        short = 0
        if os.path.isfile(os.path.join(tmp, "short.fasta")):
            short = open(os.path.join(tmp, "short.fasta")).read().count(">")
        filtered = glob.glob(os.path.join(tmp, "*filtered*"))[0]
        n_filt = open(filtered).read().count(">")
        return short, n_filt

    def test_strict_marks_deletion_short(self):
        short, n_filt = self._run("strict")
        self.assertEqual(short, 1)
        self.assertEqual(n_filt, 4)

    def test_deletion_tolerant_keeps_deletion(self):
        short, n_filt = self._run("deletion_tolerant")
        self.assertEqual(short, 0)
        self.assertEqual(n_filt, 5)

    def test_deletion_tolerant_long_uses_alignment_not_query_length(self):
        """ONT read: long query_length but acceptable alignment span -> not long."""
        region = {"chr": "ref", "start": 0, "end": 595, "name": "amplicon"}
        region_length = 595
        min_overlap = 0.95

        class FakeRead:
            reference_start = 0
            reference_end = 580
            query_alignment_length = 520
            query_length = 2500

        read = FakeRead()
        self.assertFalse(
            is_short_read(read, region, region_length, min_overlap, "deletion_tolerant")
        )
        long_default = int(region_length * (2 - min_overlap) + 200)
        self.assertFalse(is_long_read(read, long_default, "deletion_tolerant"))
        self.assertFalse(is_long_read(read, 3000, "deletion_tolerant"))
        self.assertTrue(is_long_read(read, 2000, "strict"))


if __name__ == "__main__":
    unittest.main()
