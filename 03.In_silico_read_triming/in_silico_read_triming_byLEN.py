#!/usr/bin/env python3
"""
"""
import sys, random, pysam, gzip
import numpy as np

def main():
    if len(sys.argv) < 5:
        print("Usage: samtools sort -n -o - input.bam | python3 in_silico_read_triming_byLEN.py - hist_len.txt R1.fq.gz R2.fq.gz [min_len]")
        sys.exit(1)

    bam_path = sys.argv[1]
    hist_file = sys.argv[2]
    out_r1 = sys.argv[3]
    out_r2 = sys.argv[4]
    min_len = int(sys.argv[5]) if len(sys.argv) > 5 else 50

    # 1. load hist
    with open(hist_file) as f:
        hist_tlen = [int(x.strip()) for x in f if x.strip().isdigit()]
    if not hist_tlen:
        sys.exit("no hist_len.txt file")

    # 2. process bam
    bam = pysam.AlignmentFile(bam_path if bam_path != "-" else "-", "rb")
    f1 = gzip.open(out_r1, "wt", compresslevel=6)
    f2 = gzip.open(out_r2, "wt", compresslevel=6)

    stats = {"kept": 0, "trimmed": 0, "skipped": 0}
    proc = 0
    bam_iter = iter(bam)

    try:
        for r1 in bam_iter:
            if not (r1.is_paired and r1.is_read1 and not r1.is_secondary and not r1.is_supplementary):
                continue
            try:
                r2 = next(bam_iter)
            except StopIteration:
                break
            if r2.query_name != r1.query_name:
                continue
            
            seq1, q1 = r1.get_forward_sequence(), r1.get_forward_qualities()
            seq2, q2 = r2.get_forward_sequence(), r2.get_forward_qualities()
            len1, len2 = len(seq1), len(seq2)
            
            # to trim reads
            target = random.choice(hist_tlen)
            
            if len1 <= target or len2 <= target:
                f1.write(f"@{r1.query_name} /1\n{seq1}\n+\n{''.join(chr(x+33) for x in q1)}\n")
                f2.write(f"@{r2.query_name} /2\n{seq2}\n+\n{''.join(chr(x+33) for x in q2)}\n")
                stats["kept"] += 1
            else:
                final_len = min(len1, len2, target)
                f1.write(f"@{r1.query_name} /1\n{seq1[:final_len]}\n+\n{''.join(chr(x+33) for x in q1[:final_len])}\n")
                f2.write(f"@{r2.query_name} /2\n{seq2[:final_len]}\n+\n{''.join(chr(x+33) for x in q2[:final_len])}\n")
                stats["trimmed"] += 1

            proc += 1
            if proc % 50000 == 0:
                sys.stderr.write(f"\rProcessed {proc} pairs (K:{stats['kept']} T:{stats['trimmed']})...")
                sys.stderr.flush()
    except BrokenPipeError:
        pass
    finally:
        sys.stderr.write(f"\n Done.\n")
        sys.stderr.write(f"   Total: {proc:,} pairs\n")
        sys.stderr.write(f"   Kept: {stats['kept']:,} ({stats['kept']/proc*100:.1f}%)\n")
        sys.stderr.write(f"   Trimmed: {stats['trimmed']:,} ({stats['trimmed']/proc*100:.1f}%)\n")
        f1.close(); f2.close(); bam.close()

if __name__ == "__main__":
    main()