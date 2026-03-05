#!/usr/bin/env python3
# Preprocess DNase-seq narrowPeak BED to create balanced positive/negative sets of fixed-length sequences
# Jacob Mitchell
# Usage example:
# python3 preprocess_dnase.py \
#  --bed experiment1.bed \
#  --genome /path/to/GRCh38.fa \
#  --outdir out_experiment1 \
#  --min_x 50 \
#  --min_keep_frac 0.60 \
#  --neg_ratio 1.0 \
#  --tol 0.10

import argparse
import os
import random
import statistics
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Iterable, Tuple

@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int  # half-open [start, end)

    @property
    def length(self) -> int:
        return self.end - self.start

def read_bed_3col(path: str) -> List[Interval]:
    """Read first 3 columns of a BED/narrowPeak-like file."""
    out: List[Interval] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            if end > start:
                out.append(Interval(chrom, start, end))
    return out

def sort_intervals(intervals: List[Interval]) -> List[Interval]:
    return sorted(intervals, key=lambda iv: (iv.chrom, iv.start, iv.end))

def merge_overlaps(intervals: List[Interval]) -> List[Interval]:
    """Merge overlaps/adjacent intervals per chromosome."""
    intervals = sort_intervals(intervals)
    merged: List[Interval] = []
    for iv in intervals:
        if not merged or merged[-1].chrom != iv.chrom or iv.start > merged[-1].end:
            merged.append(iv)
        else:
            last = merged[-1]
            merged[-1] = Interval(last.chrom, last.start, max(last.end, iv.end))
    return merged

def choose_x_optimize(lengths: List[int], min_keep_frac: float = 0.60,
                      min_x: int = 50, max_x: Optional[int] = None) -> int:
    """
    Choose X that is a 'middle ground':
    - keeps at least min_keep_frac of intervals (discarding short ones),
    - prefers larger X (more content captured),
    - penalizes trimming of long intervals.
    """
    if not lengths:
        raise ValueError("No lengths to optimize X.")

    L = sorted(lengths)
    n = len(L)
    if max_x is None:
        max_x = L[-1]

    def clamp(v: int) -> int:
        return max(min_x, min(int(v), max_x))

    # Candidate X values: percentiles + near-median
    cand = set()
    for p in [5,10,15,20,25,30,35,40,45,50,60,70,75,80,85,90,95]:
        idx = int(round((p/100) * (n-1)))
        cand.add(L[idx])
    med = int(statistics.median(L))
    for d in [-200,-100,-50,-25,-10,-5,0,5,10,25,50,100,200]:
        cand.add(med + d)

    candidates = sorted([c for c in cand if min_x <= c <= max_x])
    if not candidates:
        return clamp(med)

    # Precompute suffix counts and suffix sums for fast scoring:
    # for a given X, kept = count(L >= X), trim_loss = sum(L - X for L >= X) = sum_kept - kept*X
    # We'll binary-search the first index where L >= X.
    import bisect
    suffix_sum = [0]*(n+1)
    for i in range(n-1, -1, -1):
        suffix_sum[i] = suffix_sum[i+1] + L[i]

    trim_penalty = 0.25  # higher => more anti-trimming

    best_x = None
    best_score = float("-inf")

    for X in candidates:
        i = bisect.bisect_left(L, X)
        kept = n - i
        if kept / n < min_keep_frac:
            continue
        sum_kept = suffix_sum[i]
        trim_loss = sum_kept - kept * X  # total bases trimmed away
        # Score favors: many kept and larger X, penalize heavy trimming
        score = (kept * X) - (trim_penalty * trim_loss)
        if score > best_score:
            best_score = score
            best_x = X

    if best_x is None:
        # fallback: 25th percentile tends to keep many while not over-trimming
        return clamp(L[int(round(0.25*(n-1)))])
    return clamp(best_x)

def write_bed(path: str, intervals: Iterable[Interval]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for iv in intervals:
            f.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\n")

def make_positive_fixed_windows(merged: List[Interval], X: int) -> Tuple[List[Interval], int]:
    """Discard intervals shorter than X; trim longer ones to [start, start+X]."""
    pos: List[Interval] = []
    discarded = 0
    for iv in merged:
        if iv.length < X:
            discarded += 1
        else:
            pos.append(Interval(iv.chrom, iv.start, iv.start + X))
    return pos, discarded

def compute_gaps(merged: List[Interval]) -> List[Interval]:
    """Gaps between merged positive regions on each chromosome."""
    merged = sort_intervals(merged)
    gaps: List[Interval] = []
    prev_by_chr: Dict[str, Interval] = {}
    for iv in merged:
        prev = prev_by_chr.get(iv.chrom)
        if prev is not None and iv.start > prev.end:
            gaps.append(Interval(iv.chrom, prev.end, iv.start))
        prev_by_chr[iv.chrom] = iv
    return gaps

def sample_negative_windows(gaps: List[Interval], X: int, n_target: int, seed: int) -> List[Interval]:
    """
    Sample exactly-length-X windows from the *middle* of gap regions.
    No extra user parameters: hardcode a central-band strategy.
    """
    rng = random.Random(seed)
    eligible = [g for g in gaps if g.length >= X]
    if not eligible:
        return []

    # Weight by how many valid window-starts exist in the gap
    weights = [g.length - X + 1 for g in eligible]
    total_w = sum(weights)

    neg: List[Interval] = []
    for _ in range(n_target):
        # Choose a gap proportional to its available start positions
        r = rng.randrange(total_w)
        acc = 0
        chosen = eligible[-1]
        for g, w in zip(eligible, weights):
            acc += w
            if r < acc:
                chosen = g
                break

        s_min = chosen.start
        s_max = chosen.end - X  # inclusive max start
        if s_max < s_min:
            continue  # shouldn't happen due to eligibility check

        # Sample from the *central band* of [s_min, s_max]
        # Middle 50% band (hardcoded, no user-facing parameter).
        span = s_max - s_min
        if span == 0:
            start = s_min
        else:
            center = (s_min + s_max) / 2.0
            half_band = max(1, int(round(span * 0.25)))  # 25% each side => middle 50%
            lo = max(s_min, int(round(center - half_band)))
            hi = min(s_max, int(round(center + half_band)))
            if lo > hi:
                lo, hi = s_min, s_max
            start = rng.randint(lo, hi)

        neg.append(Interval(chosen.chrom, start, start + X))

    return neg

def bedtools_getfasta(genome_fa: str, bed_path: str, out_fa: str) -> None:
    subprocess.run(
        ["bedtools", "getfasta", "-fi", genome_fa, "-bed", bed_path, "-fo", out_fa],
        check=True
    )

def fasta_to_txt(fa_path: str, txt_path: str) -> int:
    """Convert FASTA to 1-seq-per-line txt (skip headers)."""
    n = 0
    with open(fa_path, "r", encoding="utf-8") as fin, open(txt_path, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            fout.write(line.upper() + "\n")
            n += 1
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True, help="Input BED/narrowPeak. Uses first 3 columns only.")
    ap.add_argument("--genome", required=True, help="Genome FASTA (e.g., GRCh38.fa)")
    ap.add_argument("--outdir", default="out", help="Output directory")

    # X optimization behavior
    ap.add_argument("--min_x", type=int, default=50, help="Minimum X allowed")
    ap.add_argument("--max_x", type=int, default=None, help="Maximum X allowed (optional)")
    ap.add_argument("--min_keep_frac", type=float, default=0.60,
                    help="Must keep at least this fraction of positives (discard short ones).")

    # Balancing
    ap.add_argument("--neg_ratio", type=float, default=1.0,
                    help="Target negatives ~= neg_ratio * positives (default 1.0)")
    ap.add_argument("--tol", type=float, default=0.10,
                    help="Allowed relative difference (e.g., 0.10 => within 10%)")

    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    raw = read_bed_3col(args.bed)
    if not raw:
        raise SystemExit("No valid intervals found in input.")

    merged = merge_overlaps(raw)
    lengths = [iv.length for iv in merged]

    X = choose_x_optimize(lengths, min_keep_frac=args.min_keep_frac,
                          min_x=args.min_x, max_x=args.max_x)

    pos, discarded = make_positive_fixed_windows(merged, X)
    n1 = len(pos)
    if n1 == 0:
        raise SystemExit(f"Chosen X={X} discards all positives. Try lowering --min_x or --min_keep_frac.")

    # negatives: aim close (within tol). Default is roughly equal.
    n2_target = max(1, int(round(args.neg_ratio * n1)))
    gaps = compute_gaps(merged)
    neg = sample_negative_windows(gaps, X, n2_target, seed=args.seed)

    # If we couldn't sample enough negatives (rare, but possible), we just take what we can.
    n2 = len(neg)
    rel_diff = abs(n2 - n1) / n1

    # Write BEDs
    pos_bed = os.path.join(args.outdir, "positive.bed")
    neg_bed = os.path.join(args.outdir, "negative.bed")
    write_bed(pos_bed, pos)
    write_bed(neg_bed, neg)

    print(f"Chosen X = {X}")
    print(f"Positives kept: {n1} (discarded short: {discarded})")
    print(f"Negatives sampled: {n2} (target: {n2_target}, rel_diff_vs_pos: {rel_diff:.3f})")
    if rel_diff > args.tol:
        print(f"WARNING: negative/positive counts differ by > {args.tol*100:.0f}%. "
              f"Likely not enough eligible gap space at length X.")

    # Produce sequences
    pos_fa = os.path.join(args.outdir, "positive.fa")
    neg_fa = os.path.join(args.outdir, "negative.fa")
    bedtools_getfasta(args.genome, pos_bed, pos_fa)
    bedtools_getfasta(args.genome, neg_bed, neg_fa)

    pos_txt = os.path.join(args.outdir, "positive.txt")
    neg_txt = os.path.join(args.outdir, "negative.txt")
    npos = fasta_to_txt(pos_fa, pos_txt)
    nneg = fasta_to_txt(neg_fa, neg_txt)

    print(f"Wrote {pos_txt} ({npos} sequences)")
    print(f"Wrote {neg_txt} ({nneg} sequences)")

if __name__ == "__main__":
    main()