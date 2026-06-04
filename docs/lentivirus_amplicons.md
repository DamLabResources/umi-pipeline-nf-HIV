# Lentivirus / SIVmac239m2 UMI amplicons

This fork adds optional behavior for integrated lentivirus amplicons with large internal deletions.

## Reference and BED

Use a provirus reference (e.g. `sivmac239m2.ltr_masked.fa`) and one or more BED targets:

```text
SIVmac239m2	1179	1774	amplicon
SIVmac239m2	251	6900	amplicon
```

Short and long targets can be combined in one BED file for multi-target runs.

## Split-read filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_overlap` | `0.95` | Fraction of BED span required for retention |
| `split_read_filter_mode` | `strict` | `strict`: short = `query_alignment_length`; long = full `query_length` vs BED span. `deletion_tolerant`: short = **reference overlap** with BED; long = **`query_alignment_length`** (not full read length), so ONT reads with internal deletions are not mis-binned as `long` |

Enable deletion tolerance when deletion-rich reads are classified as `short` in verbose stats (`reads_short`).

Set `max_query_length` high enough for Nanopore read lengths (e.g. **3000** for ~595 bp targets, **10000** for ~6.5 kb SLN targets). Without it, the default long threshold (~`region_length + 5%` + adapters) bins almost all ONT reads as `long` even in `deletion_tolerant` mode.

## Length-based subclustering (same UMI, different haplotype length)

After UMI-tag subclustering (`max_dist_umi`), optionally split again by full amplicon length before Medaka polishing:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cluster_by_amplicon_length` | `false` | Enable length splitting |
| `max_amplicon_length_diff_bp` | `null` | Max absolute length difference within one polishing group |
| `max_amplicon_length_diff_pct` | `null` | Max percent length difference within one polishing group |

Two reads stay in the same polishing group if **either** threshold is satisfied (bp **or** %). Set at least one threshold when the feature is enabled.
Biological deletions are typically large, nanopore sequencing can create small deletions.

Suggested starting points for length-based subclustering:

- `max_amplicon_length_diff_bp = 80–150`
- `max_amplicon_length_diff_pct = 8–12`


## Other recommended settings

- `reference_based_polishing = false` (POA) for deletion-heavy data
- `balance_strands = false` if reverse-strand clusters are sparse
- `min_length = 16` when using 17–19 bp UMI designs
- `use_context = true` with primer sequences flanking UMIs when tags are repetitive

## Example config

See [`config/examples/sivmac239m2.config`](../config/examples/sivmac239m2.config).

## Smoke test

```bash
nextflow run . -profile test,docker -config config/sivmac_smoke.config
```

Fixtures are generated from SIVmac239m2 with:

```bash
python3 scripts/build_sivmac_test_fixtures.py
```

