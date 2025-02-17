#!/bin/env python

import os
import sys
import pysam
import logging
import argparse
import matplotlib

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib import ticker
from typing import Generator

matplotlib.use("Agg")

logging.basicConfig(
    format="%(levelname)s (%(asctime)s): %(message)s (Line: %(lineno)d [%(filename)s])",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.WARNING,
)

OUTPUT_HEADER = ["Sample", "Coverage(X)", "Coverage(X,100kb+)", "N50(kb)", "Reads"]


def NX(nums, X):
    nums = sorted(nums, key=int, reverse=True)
    datathresh = sum(nums) * (X / 100.0)
    total = 0
    for num in nums:
        total += num
        if total >= datathresh:
            return num
    return 0


def get_n50(vals: np.ndarray, *, sort: bool = False):
    if sort:
        vals = np.flip(np.sort(vals))
    vals_csum = np.cumsum(vals)
    val_mdpt = vals_csum[-1] // 2
    return vals[np.sum(vals_csum <= val_mdpt)] / 1000


def gb_formatter(x, pos):
    return f"{x / 1000000000} Gbp"


def mb_formatter(x, pos):
    return f"{int(x / 1000000)} Mbp"


def kb_formatter(x, pos):
    return f"{int(x / 1000)} kbp"


def load_sum_hist(len_list=None, bins=None):
    bin_df = pd.DataFrame({"BIN": bins})
    bin_df["BIN"] = bin_df.apply(lambda x: x // BIN_SIZE)
    len_df = pd.DataFrame({"SUM": len_list})
    len_df["BIN"] = len_df["SUM"].apply(lambda x: x // BIN_SIZE)
    sum_len_df = len_df.groupby("BIN").sum()
    sum_hist = pd.merge(bin_df, sum_len_df, on="BIN", how="left").fillna(0)
    return sum_hist


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    "bams",
    nargs="+",
    required=True,
    help="Unmapped BAM files.",
)
parser.add_argument("--sm", "-s", default=None, help="Sample name.")
parser.add_argument(
    "--genome", "-g", type=float, required=False, help="Genome size in Gbp", default=3.1
)
parser.add_argument(
    "--outfile",
    "-o",
    type=argparse.FileType("at"),
    required=False,
    help="Output file to write to",
    default=sys.stdout,
)
parser.add_argument(
    "--tab",
    "-t",
    required=False,
    action="store_true",
    help="Output information in a tab-delimited format",
)
parser.add_argument(
    "--plot",
    "-p",
    type=str,
    required=False,
    help="Plot the read length distribution and save to the provided argument",
)
parser.add_argument(
    "--log",
    "-l",
    required=False,
    action="store_true",
    default=False,
    help="Plot the read length distribution on a log scale",
)
parser.add_argument(
    "--length_limit",
    "-x",
    type=int,
    required=False,
    help="Limit of the read length to be displayed on the plot",
)
parser.add_argument(
    "--cumulative",
    "-u",
    required=False,
    action="store_true",
    default=False,
    help="Plot the cumulative sum of base pairs by read length",
)
parser.add_argument(
    "--window",
    "-w",
    required=False,
    type=int,
    default=10000,
    help="Read length window size for plotting bins",
)
parser.add_argument(
    "--plot_type", choices=["original", "cdf"], default="cdf", help="Plot type."
)
args = parser.parse_args()


def read_bam(infile: str) -> Generator[tuple[str, int], None, None]:
    with pysam.AlignmentFile(infile, check_sq=False) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or (not rec.is_secondary and not rec.is_supplementary):
                yield rec.query_name, len(rec.query_sequence)


lengths = np.array(length for bam in args.bams for _, length in read_bam(bam))

if args.tab:
    args.outfile.write("\t".join(OUTPUT_HEADER) + "\n")

len_list = np.flip(np.sort(lengths))
len_list_k = lengths[lengths >= 100_000]
n50 = get_n50(len_list)

coverage = np.sum(len_list) / (args.genome * 1_000_000_000)
coverage_k = np.sum(len_list_k) / (args.genome * 1_000_000_000)

sm = args.sample
if not args.tab:
    args.outfile.write(f"Sample: {sm}\n")
    args.outfile.write(
        "Coverage (X): {:,.2f}\n"
        "Coverage 100k+ (X): {:,.2f}\n"
        "N50 (kbp):   {:,.2f}\n"
        "Reads:	 {:,.0f}\n".format(
            coverage,
            coverage_k,
            n50,
            len(len_list),
        )
    )
else:
    out_df = pd.DataFrame.from_dict(
        dict(
            zip(
                OUTPUT_HEADER,
                [
                    [sm],
                    [round(coverage, 2)],
                    [round(coverage_k, 2)],
                    [round(n50, 2)],
                    [len(len_list)],
                ],
            )
        )
    )
    out_df.to_csv(args.outfile, sep="\t", index=False, header=False, mode="a")

if args.plot_dir:
    os.makedirs(args.plot_dir, exist_ok=True)

    if args.plot_type == "cdf":
        sns.set(font_scale=3)
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(20, 12))

        log_lengths = np.log10(np.clip(lengths, 10, None))

        # make histogram
        sns.distplot(log_lengths, bins=100, kde=False, rug=False, ax=ax)

        # get maxes
        myy = plt.gca().get_ylim()[1]
        mymax = max(lengths)

        # add vertical lines
        vals = [
            np.median(lengths),
            lengths.mean(),
            NX(lengths, 50.0),
            NX(lengths, 1.0),
            mymax,
        ]
        names = ["Median", "Mean", "N50", "N1", "Max"]
        divs = [0.9, 0.8, 0.7, 0.6, 0.5]
        for name, val, div in zip(names, vals, divs):
            plt.axvline(x=np.log10(val), color="darkred", linestyle="dashed")
            label = "{}={}".format(name, round(val / 1000, 1))
            plt.text(np.log10(val) + 0.01, myy * div, label, fontsize=28)

        xts = [2, 3, 4, 5, 6, 7]
        xls = ["0.1", "1", "10", "100", "1,000", "10,000"]

        # plot
        plt.xticks(xts, xls)
        Gb = sum(lengths) / 1000000000
        plt.xlabel("Length (kbp), Total Gb = {:.2f}".format(Gb))
        plt.ylabel("Number of reads")
        plt.tight_layout()
        plt.savefig(args.plot, bbox_inches="tight")
    else:
        BIN_SIZE = args.window
        fig, ax = plt.subplots()

        if args.log:
            ax.set_xscale("log")
            ax.set_xlabel("Read Length (log)")
        else:
            ax.set_xlabel("Read Length")

        bins = np.arange(0, max(len_list) + BIN_SIZE, BIN_SIZE)
        sum_hist = load_sum_hist(len_list=len_list, bins=bins)[:-1]

        ax.xaxis.set_major_formatter(ticker.FuncFormatter(kb_formatter))

        if args.cumulative:
            if args.length_limit:
                ax.set_xlim(args.length_limit, 0)
            else:
                ax.set_xlim(max(bins), 0)

            cumulative_hist = np.cumsum(sum_hist[::-1])[::-1]
            ax.bar(bins[:-1], cumulative_hist["SUM"], width=BIN_SIZE, align="edge")
            ax.set_ylabel("Cumulative Sum of the Number of Bases")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(gb_formatter))
        else:
            if args.length_limit:
                ax.set_xlim(0, args.length_limit)

            ax.bar(bins[:-1], sum_hist["SUM"], width=BIN_SIZE, align="edge")
            ax.set_ylabel("Sum of the Number of Bases")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(mb_formatter))

        plt.title(sm)
        plt.savefig(args.plot, bbox_inches="tight")
