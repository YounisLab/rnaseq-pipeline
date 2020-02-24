#!/usr/bin/env python

"""
Remove transgene entries from a regtools junctions.bed file, given
a reference .bed file containing gene locus.
"""

import sys
import argparse
import os.path
import subprocess as sp

parser = argparse.ArgumentParser(
    description="Remove transgene entries from a regtools junctions.bed file, given \
    a reference .bed file containing gene locus.",
    epilog="Usage eg: python remove_transgene.py <hg38.refGene_gene_longest.bed> <junctions.bed> <output.bed>"
)

parser.add_argument("refgene_path", help="Path to reference .bed file containing gene loci.")
parser.add_argument("junctions_path", help="Path to .bed file containing junctions from regtools.")
parser.add_argument("outname", help="Full path to output file.")

args = parser.parse_args()

if not os.path.exists(args.refgene_path):
    print("Reference .bed file does not exist! Exiting...")
    exit(1)
if not os.path.exists(args.junctions_path):
    print(".bed file containing junctions does not exist! Exiting...")
    exit(1)

refgene_bed = open(args.refgene_path)

print(">> Running intersectBed to determine junction spans.")
cmd = "intersectBed -a %s -b %s -wao -f 1" % (args.junctions_path, args.refgene_path)
intersect_f = open(args.outname + "_intersect.bed", "w+")
sp.call(cmd.split(" "), stdout=intersect_f)

print(">> Filtering transgene junctions...")
intersect_f.seek(0)
lines = intersect_f.readlines()
total = len(lines)
transgene_count = 0
unknown_strand_count = 0
tmp_out_f = open(args.outname + "_tmp", "w+")
for progress, line in enumerate(lines):
#    import pdb; pdb.set_trace()
    sys.stdout.write( "\r>> Row processing progress: %f percent completed..." % (float(progress + 1)/total * 100))
    sys.stdout.flush()
    splitted = line.split("\t")
    if splitted[-1].strip() == "0":
        transgene_count += 1
        continue
    if splitted[5].strip() == "?":
        unknown_strand_count += 1
        splitted[5] = "+"
    
    tmp_out_f.write("\t".join(splitted[:12]) + "\n")

tmp_out_f.close()
sys.stdout.write("\n")
print(">> Filtering complete. Transgene count: %d | Unknown strand count: %d" % \
    (transgene_count, unknown_strand_count))

# Since some junctions can map to multiple genes,
# we need to uniquely extract rows.
print(">> Removing duplicates...")
cmd = "awk !seen[$0]++ %s" % (args.outname + "_tmp")
out_f = open(args.outname, "w+")
sp.call(cmd.split(" "), stdout=out_f)
out_f.close()
print(">> Cleaning up and exiting gracefully...")
os.remove(args.outname + "_tmp")
exit(0)
