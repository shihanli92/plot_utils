#!/usr/bin/env python3
"""Download the latest VDJdb release and extract fully-paired human TCR records.

Resolves the latest release of antigenomics/vdjdb-db via the GitHub API (no
hardcoded tag), downloads the data zip attached to the release, and extracts the
paired-chain table (vdjdb_full.txt). Keeps only fully-paired human records --
rows where both cdr3.alpha and cdr3.beta are populated -- and writes
vdjdb_paired.csv.

If no *_full.txt is present, falls back to the slim table and reconstructs pairs
by grouping on complex.id != 0.
"""
import io
import zipfile

import pandas as pd
import requests

GITHUB_API = "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest"
HUMAN = "HomoSapiens"
OUTPUT = "vdjdb_paired.csv"


def resolve_release():
    """Return (tag, zip_download_url) for the latest release."""
    rel = requests.get(GITHUB_API, timeout=60).json()
    zip_assets = [a for a in rel["assets"] if a["name"].lower().endswith(".zip")]
    if not zip_assets:
        raise RuntimeError("No zip asset found on the latest VDJdb release")
    asset = zip_assets[0]
    return rel["tag_name"], asset["name"], asset["browser_download_url"]


def is_human(df):
    """Boolean mask selecting human TCR records (species column)."""
    return df["species"].astype(str).str.strip() == HUMAN


def from_full(zf, member):
    """Extract fully-paired human records from the full (paired) table."""
    df = pd.read_csv(zf.open(member), sep="\t", dtype=str, low_memory=False)
    df = df[is_human(df)]
    paired = (
        df["cdr3.alpha"].notna() & (df["cdr3.alpha"].str.strip() != "")
        & df["cdr3.beta"].notna() & (df["cdr3.beta"].str.strip() != "")
    )
    return df[paired].reset_index(drop=True)


def from_slim(zf, member):
    """Reconstruct paired human records from the slim table via complex.id."""
    df = pd.read_csv(zf.open(member), sep="\t", dtype=str, low_memory=False)
    df = df[is_human(df)]
    # complex.id == 0 means an unpaired/single-chain record; drop those.
    df = df[df["complex.id"].fillna("0").str.strip() != "0"]
    rows = []
    gene_col = "gene"  # slim table stores chain identity (TRA/TRB) here
    for cid, grp in df.groupby("complex.id"):
        alpha = grp[grp[gene_col] == "TRA"]
        beta = grp[grp[gene_col] == "TRB"]
        if alpha.empty or beta.empty:
            continue  # not a complete pair
        a = alpha.iloc[0]
        b = beta.iloc[0]
        shared = a.drop(labels=["gene", "cdr3", "v.segm", "j.segm"], errors="ignore").to_dict()
        shared.update({
            "complex.id": cid,
            "cdr3.alpha": a["cdr3"], "v.alpha": a.get("v.segm"), "j.alpha": a.get("j.segm"),
            "cdr3.beta": b["cdr3"], "v.beta": b.get("v.segm"), "j.beta": b.get("j.segm"),
        })
        rows.append(shared)
    return pd.DataFrame(rows).reset_index(drop=True)


def main():
    tag, name, url = resolve_release()
    print(f"Latest VDJdb release: {tag} -> {name}")
    data = requests.get(url, timeout=600).content
    zf = zipfile.ZipFile(io.BytesIO(data))
    full = [n for n in zf.namelist() if n.endswith("_full.txt")]
    slim = [n for n in zf.namelist() if n.endswith(".slim.txt")]
    if full:
        print(f"Using paired-chain table: {full[0]}")
        result = from_full(zf, full[0])
    elif slim:
        print(f"No *_full.txt found; reconstructing pairs from slim table: {slim[0]}")
        result = from_slim(zf, slim[0])
    else:
        raise RuntimeError("Neither a *_full.txt nor a *.slim.txt table was found")
    result.to_csv(OUTPUT, index=False)
    print(f"\nWrote {OUTPUT}")
    print(f"Final row count: {len(result)}")
    print(f"Columns ({len(result.columns)}): {list(result.columns)}")


if __name__ == "__main__":
    main()
