#!/usr/bin/env python
"""Compare two output trees CSV-by-CSV and report numeric agreement.

Usage:
  python code/compare_outputs.py --a <shipped-root> --b <rerun-root> [--tol 1e-4]

Reports one row per file. A file whose shared numeric columns are empty is
flagged NO_NUMERIC_OVERLAP rather than silently passing.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--a', required=True, help='reference tree (e.g. the package root)')
ap.add_argument('--b', required=True, help='rerun tree (e.g. py_outputs)')
ap.add_argument('--tol', type=float, default=1e-4)
ap.add_argument('--out', default=None, help='optional CSV path for the report')
args = ap.parse_args()

A, B = Path(args.a).resolve(), Path(args.b).resolve()
rows, worst = [], 0.0

subdirs = ['height/outputs', 'bmi/outputs', 'code']
for sub in subdirs:
    src = A / sub
    if not src.is_dir():
        continue
    for f in sorted(src.glob('*.csv')):
        rel = f'{sub}/{f.name}'
        other = B / sub / f.name
        if not other.exists():
            rows.append(dict(file=rel, status='MISSING_IN_B', n=0, max_abs=np.nan, max_rel=np.nan, note=''))
            continue
        try:
            da, db = pd.read_csv(f), pd.read_csv(other)
        except Exception as e:
            rows.append(dict(file=rel, status='READ_ERROR', n=0, max_abs=np.nan, max_rel=np.nan, note=str(e)[:60]))
            continue
        note = ''
        if len(da) != len(db):
            rows.append(dict(file=rel, status='ROW_COUNT', n=0, max_abs=np.nan, max_rel=np.nan,
                             note=f'{len(da)} vs {len(db)}'))
            continue
        shared = [c for c in da.columns if c in db.columns]
        only_a = [c for c in da.columns if c not in db.columns]
        if only_a:
            note = 'cols_missing_in_B=' + ';'.join(only_a)
        n, max_abs, max_rel = 0, 0.0, 0.0
        text_diff = False
        for c in shared:
            x, y = da[c], db[c]
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                xv, yv = x.to_numpy(float), y.to_numpy(float)
                ok = np.isfinite(xv) & np.isfinite(yv)
                if ok.any():
                    d = np.abs(xv[ok] - yv[ok])
                    denom = np.maximum(np.maximum(np.abs(xv[ok]), np.abs(yv[ok])), 1e-12)
                    n += int(ok.sum()); max_abs = max(max_abs, float(d.max()))
                    max_rel = max(max_rel, float((d / denom).max()))
            elif not x.astype(str).equals(y.astype(str)):
                text_diff = True
        if n == 0:
            status = 'NO_NUMERIC_OVERLAP'
        elif max_abs > args.tol:
            status = 'NUMERIC_DIFFERS'
        elif text_diff:
            status = 'TEXT_DIFFERS'
        else:
            status = 'OK'
        worst = max(worst, max_abs)
        rows.append(dict(file=rel, status=status, n=n, max_abs=max_abs, max_rel=max_rel, note=note))

rep = pd.DataFrame(rows)
pd.set_option('display.width', 200, 'display.max_colwidth', 60)
print(rep.to_string(index=False))
bad = rep[~rep.status.isin(['OK'])]
print(f'\nfiles compared: {len(rep)}   clean: {(rep.status == "OK").sum()}   flagged: {len(bad)}')
print(f'largest absolute difference anywhere: {worst:.3e}  (tolerance {args.tol:g})')
if args.out:
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.out, index=False)
    print(f'wrote {args.out}')
sys.exit(0)
