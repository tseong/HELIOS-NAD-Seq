#!/usr/bin/env python3
import sys
from pathlib import Path
import re
from itertools import groupby

"""
Usage:
    python make_intergenic_gtf.py input.gtf intergenic.gtf

Emits one intergenic feature per gene, preserving original GTF fields.
"""

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

in_gtf, out_gtf = Path(sys.argv[1]), Path(sys.argv[2])

# 0) slurp and write headers
header_lines = []
with in_gtf.open() as fin:
    for line in fin:
        if line.startswith('#'):
            header_lines.append(line)
        else:
            break

# Regex to parse GTF attributes
_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
def parse_attrs(s):
    return {k: v for k, v in _attr_re.findall(s)}

# 1) collect all gene records (keeping all original columns we?~@~Yll need)
genes = []
with in_gtf.open() as fin:
    for L in fin:
        if L.startswith('#'):
            continue
        cols = L.rstrip('\n').split('\t')
        if len(cols) < 9 or cols[2] != 'gene':
            continue
        chrom, s, e, src, score, strand, frame, attr = (
            cols[0],
            int(cols[3]),
            int(cols[4]),
            cols[1],
            cols[5],
            cols[6],
            cols[7],
            cols[8],
        )
        genes.append((chrom, s, e, src, score, strand, frame, attr))

# sort and build a map of gene_id -> next_gene_start on each chrom
genes.sort(key=lambda x: (x[0], x[1]))
_next_start = {}
for chrom, grp in groupby(genes, key=lambda x: x[0]):
    lst = list(grp)
    for i, (c, s, e, src, score, strand, frame, attr) in enumerate(lst):
        gid = parse_attrs(attr).get('gene_id')
        _next_start[gid] = lst[i+1][1] if i < len(lst)-1 else None
      # 2) write out the intergenic GTF
with out_gtf.open('w') as fout:
    # write original headers
    for hl in header_lines:
        fout.write(hl)

    prev_end = {}  # track last gene end per chrom
    for chrom, start, end, src, score, strand, frame, attr in genes:
        gid = parse_attrs(attr).get('gene_id','.')

        if strand == '+':
            ig_end   = end
            pe       = prev_end.get(chrom, 0)
            ig_start = max(pe+1, start-100)
            # if the ?~@~\previous?~@~] gene actually overlaps this one:
            if pe >= start:
                ig_start = start

        else:  # strand == '-'
            ig_start = start
            ns       = _next_start.get(gid)
            if ns is None:
                ig_end = end + 100
            else:
                gap = ns - end - 1
                if gap <= 0:
                    # overlapping downstream gene
                    ig_end = end
                elif gap <= 100:
                    ig_end = ns - 1
                else:
                    ig_end = end + 100

        # emit the intergenic line, preserving cols 1,2,6,7,9
        out_cols = [
            chrom,
            src,
            'intergenic',
            str(ig_start),
            str(ig_end),
            score,
            strand,
            frame,
            attr
        ]
        fout.write('\t'.join(out_cols) + '\n')
        prev_end[chrom] = end
