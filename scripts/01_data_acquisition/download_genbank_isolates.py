#!/usr/bin/env python3
import csv, os, re, subprocess, time
from collections import defaultdict

PROJECT = os.path.expanduser("~/kp_liver_project")
infile = os.path.join(PROJECT, "lists/genbank_map.csv")
outdir = os.path.join(PROJECT, "genomes_fasta/genbank_isolates")
logdir = os.path.join(PROJECT, "logs")
os.makedirs(outdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)

def safe_name(s: str) -> str:
    s = str(s).strip()
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)

ok_log = os.path.join(logdir, "genbank_download_ok.tsv")
fail_log = os.path.join(logdir, "genbank_download_fail.tsv")

# Fix Excel BOM in header if present
with open(infile, "rb") as f:
    raw = f.read()
raw = raw.lstrip(b"\xef\xbb\xbf")
with open(infile, "wb") as f:
    f.write(raw)

rows = []
with open(infile, newline="") as f:
    reader = csv.DictReader(f)
    # expected column names
    gid_key = "Genome ID"
    acc_key = "GenBank Accessions"
    if gid_key not in reader.fieldnames or acc_key not in reader.fieldnames:
        raise SystemExit(f"Unexpected columns: {reader.fieldnames}")

    for r in reader:
        gid = (r[gid_key] or "").strip()
        accs = (r[acc_key] or "").strip().strip('"')
        if not gid or not accs:
            continue
        acc_list = [a.strip() for a in accs.split(",") if a.strip()]
        if acc_list:
            rows.append((gid, acc_list))

gid_seen = defaultdict(int)

with open(ok_log, "w") as ok, open(fail_log, "w") as fail:
    ok.write("GenomeID\tAccessions\tOutputFile\n")
    fail.write("GenomeID\tAccession\tReason\n")

    for gid, acc_list in rows:
        gid_seen[gid] += 1
        tag = gid if gid_seen[gid] == 1 else f"{gid}_dup{gid_seen[gid]}"
        outpath = os.path.join(outdir, safe_name(tag) + ".fasta")

        # overwrite output for repeatability
        open(outpath, "w").close()

        any_success = False
        for acc in acc_list:
            cmd = ["efetch", "-db", "nuccore", "-id", acc, "-format", "fasta"]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0 or not res.stdout.strip():
                fail.write(f"{gid}\t{acc}\tEfetch_failed\n")
            else:
                with open(outpath, "a") as out:
                    out.write(res.stdout)
                    if not res.stdout.endswith("\n"):
                        out.write("\n")
                any_success = True
            time.sleep(0.34)  # be nice to NCBI

        if any_success:
            ok.write(f"{gid}\t{','.join(acc_list)}\t{outpath}\n")
        else:
            os.remove(outpath)

print(f"Wrote GenBank isolate FASTAs to: {outdir}")
print(f"OK log:   {ok_log}")
print(f"FAIL log: {fail_log}")
