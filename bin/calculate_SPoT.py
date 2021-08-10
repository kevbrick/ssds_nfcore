#!/usr/bin/env python
from collections import defaultdict
from pybedtools import BedTool
import argparse
import os
import sys
import math

parser = argparse.ArgumentParser(
    description="Calculate the signal portion of tags / fragments (SPoT)."
)

## Allow "Just show version" mode
parser.add_argument("--reads_bed", help="Reads/Fragments BED")
parser.add_argument("--reads_bam", help="Reads BAM")
parser.add_argument("--intervals_bed", help="Genomic intervals BED")
parser.add_argument("--name", help="Reads name")
parser.add_argument("--iname", help="Intervals name")
parser.add_argument("--g", help="Genome index file")
parser.add_argument("--o", default="no_name", help="Output prefix")
parser.add_argument(
    "--rand", action="store_true", help="Also show SPoT for a randomized set"
)
parser.add_argument("--vers", action="store_true", help="Echo version")
parser.add_argument("--v", action="store_true", help="Verbose mode")

args = parser.parse_args()

if args.vers == True:
    print("1.0.0")
    quit()

if args.reads_bed is not None:
    reads = BedTool(args.reads_bed)
    intervals_bed = BedTool(args.intervals_bed).sort()
    rand_intervals_bed = (
        BedTool(args.intervals_bed).shuffle(g=args.g, chrom=True, seed=42).sort()
    )
elif args.reads_bam is not None:
    reads = BedTool(args.reads_bam)
    intervals_bed = BedTool(args.intervals_bed).sort(g=args.g)
    rand_intervals_bed = (
        BedTool(args.intervals_bed)
        .shuffle(g=args.g, chrom=True, seed=42)
        .sort(g=args.g)
    )
else:
    sys.exit("Reads BED / BAM not provided")

out_report = open(args.o + ".SSDS_SPoT_report.txt", "w")

at_intervals = reads.intersect(b=intervals_bed, u=True, sorted=True)
interval_SPoT = round(float(at_intervals.count()) / float(reads.count()), 3)
out_report.write(args.name + "_SPoT\t" + args.iname + "\t" + str(interval_SPoT) + "\n")

if args.rand:
    at_random = reads.intersect(b=rand_intervals_bed, u=True, sorted=True)
    random_SPoT = round(float(at_random.count()) / float(reads.count()), 3)
    out_report.write(
        args.name + "_SPoT\t" + args.iname + "(R)" + "\t" + str(random_SPoT) + "\n"
    )
