"""Run TIDE locally with tidepy (official Python implementation by the TIDE authors).

Usage: python tide_run.py <input.txt> <output.csv>
input.txt: tab-separated, genes (symbols) x samples, log2(TPM+1) already centred by gene (row mean subtracted).
"""
import sys

import pandas as pd
from tidepy.pred import TIDE

expr = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
res = TIDE(expr, cancer="Other", pretreat=False, vthres=0.0, ignore_norm=True)
res.index.name = "Patient"
res.to_csv(sys.argv[2])
print(f"TIDE done: {res.shape[0]} samples")
