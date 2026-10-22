#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:45

Open the all_accessions list from efetch_to_tsv.py and get the taxonomic
information for every assembly accession.

The taxonomy is resolved with only TWO batched Entrez queries, no matter how
many accessions are requested:

  1. efetch -db assembly -input <chunk> -format docsum
         assembly accession -> taxid (element <Taxid> of the docsum)
  2. efetch -db taxonomy -input <chunk> -format xml
         taxid -> lineage (the <LineageEx> of each record)

The previous implementation ran, *per accession*, a pipeline of
`esearch | efetch | elink | esearch | efetch | efetch` (six network requests)
with a fixed sleep between accessions. On a run with 3,511 accessions that
took ~4 h 19 min (with ~21k transient `curl (56)` SSL errors, each retried by
edirect). The batched version resolves the same input in well under a minute
because the number of requests drops from 6*N to 2*ceil(N/batch).

Failures are handled per chunk with retries and recursive splitting, so a
transient NCBI error only affects a handful of accessions (which are written
as "NA" rows), never the whole run.

Output format (unchanged, consumed by taxonomy_itol_files.py and
gene_phylo_distribution.py):

    assembly    taxid    scientific_name    rank

One row per ancestor in the lineage of the assembly taxid. Names and ranks are
reported verbatim by NCBI (spaces replaced by "_"); nothing is renamed or
merged here.

@author: mlonigro
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

# Number of identifiers sent to a single efetch call. Entrez accepts long
# `-input` files through POST, so these are limited only by the response size
# (a few MB per chunk), not by the command-line length.
ASSEMBLY_BATCH = 400
TAXONOMY_BATCH = 300

# Retries for a failing chunk before it is split in half. edirect already
# retries internally, so a handful of outer attempts is enough to ride out the
# occasional `curl (56) OpenSSL SSL_read: unexpected eof` responses.
MAX_ATTEMPTS = 3
RETRY_SLEEP_SECONDS = (5, 10)

# A chunk is never split below this size: at that point the identifiers are
# reported as failed instead of hammering NCBI with tiny requests.
MIN_SPLIT_SIZE = 10

HEADER = "assembly\ttaxid\tscientific_name\trank\n"

# Any accession mentioned inside an assembly docsum (current accession, last
# major release, GCA<->GCF synonym and FTP paths) belongs to that assembly.
ACCESSION_RE = re.compile(r"GC[AF]_\d+(?:\.\d+)?")


def base_accession(accession):
    """GCA_000143165.1 -> GCA_000143165, so version bumps still match."""
    return accession.rsplit(".", 1)[0]


def parse_xml(text):
    """Parse an efetch response, tolerating concatenated documents.

    NCBI occasionally returns more than one XML document in a single response
    (the XML declaration and DOCTYPE are then repeated in the middle of the
    body), which makes ElementTree fail with "junk after document element".
    Wrapping the response in a synthetic root recovers every document.
    """
    try:
        return ET.fromstring(text)
    except ET.ParseError:
        cleaned = re.sub(r"<\?xml[^>]*\?>", "", text)
        cleaned = re.sub(r"<!DOCTYPE[^>]*>", "", cleaned)
        return ET.fromstring(f"<syntenitol_root>{cleaned}</syntenitol_root>")


def sanitize(value):
    """NCBI names/ranks are written as-is, only spaces become underscores."""
    value = (value or "").strip()
    return value.replace(" ", "_") if value else "NA"


def efetch_xml(db, ids, fmt):
    """Run one `efetch` call and return the parsed XML root.

    Identifiers are passed through a temporary file (`-input`) so the command
    line never grows with the number of queries. Any failure (non-zero exit,
    empty body, invalid XML) raises RuntimeError so the caller can retry/split.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tmp:
        tmp.write("\n".join(ids) + "\n")
        input_path = tmp.name
    try:
        proc = subprocess.run(
            ["efetch", "-db", db, "-input", input_path, "-format", fmt],
            capture_output=True,
            text=True,
        )
    finally:
        os.remove(input_path)

    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip().splitlines()
        detail = stderr[-1] if stderr else "no stderr"
        raise RuntimeError(f"efetch -db {db} exited with {proc.returncode}: {detail}")
    if not proc.stdout.strip():
        raise RuntimeError(f"efetch -db {db} returned an empty response")
    try:
        return parse_xml(proc.stdout)
    except ET.ParseError as exc:
        raise RuntimeError(f"efetch -db {db} returned invalid XML: {exc}")


def parse_assembly_docsum(root):
    """DocumentSummarySet XML -> {assembly_accession: taxid}.

    Both the current accession and the GenBank/RefSeq synonyms are indexed,
    because NCBI answers a query for a GCA accession with the paired GCF record
    (same assembly, same taxid). Accessions are also indexed without their
    version suffix, so a query for an older version still resolves.
    """
    mapping = {}
    for summary in root.iter("DocumentSummary"):
        taxid = (summary.findtext("Taxid") or "").strip()
        if not taxid:
            continue
        # Serializing the record and extracting every accession it mentions
        # covers <AssemblyAccession>, <LastMajorReleaseAccession>,
        # <Synonym><Genbank>/<RefSeq> and the <FtpPath_*> entries.
        serialized = ET.tostring(summary, encoding="unicode")
        accessions = set(ACCESSION_RE.findall(serialized))
        current = (summary.findtext("AssemblyAccession") or "").strip()
        if current:
            accessions.add(current)
        for accession in accessions:
            mapping.setdefault(accession, taxid)
            mapping.setdefault(base_accession(accession), taxid)
    return mapping


def parse_taxonomy_xml(root):
    """TaxaSet XML -> {taxid: [(lineage_taxid, name, rank), ...]}.

    Only the <LineageEx> is kept (the ancestors of the requested taxid), which
    is exactly what the previous implementation wrote: one row per ancestor.
    """
    if root.tag == "Taxon":
        taxa = [root]
    else:
        # <TaxaSet><Taxon>... (normal case) or a synthetic root with several
        # concatenated documents.
        taxa = root.findall("Taxon")
        if not taxa:
            for taxa_set in root.iter("TaxaSet"):
                taxa.extend(taxa_set.findall("Taxon"))

    mapping = {}
    for taxon in taxa:
        taxid = (taxon.findtext("TaxId") or "").strip()
        if not taxid:
            continue
        lineage = []
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is not None:
            for node in lineage_ex.findall("Taxon"):
                lineage.append((
                    (node.findtext("TaxId") or "").strip() or "NA",
                    sanitize(node.findtext("ScientificName")),
                    sanitize(node.findtext("Rank")),
                ))
        mapping[taxid] = lineage
    return mapping


def fetch_and_parse(ids, db, fmt, parser, label):
    """Fetch `ids` in chunks and parse every response.

    Returns (parsed_mapping, failed_ids). A failing chunk is retried with
    exponential backoff; if it still fails it is split in half and each half is
    retried independently, until MIN_SPLIT_SIZE. This keeps a single bad
    response from discarding hundreds of accessions.
    """
    parsed = {}
    failed = []
    total = len(ids)
    batch_size = ASSEMBLY_BATCH if db == "assembly" else TAXONOMY_BATCH
    n_batches = (total + batch_size - 1) // batch_size

    def process(chunk, batch_label):
        for attempt in range(MAX_ATTEMPTS):
            try:
                root = efetch_xml(db, chunk, fmt)
                parsed.update(parser(root))
                return
            except Exception as exc:  # noqa: BLE001 - any NCBI/XML failure is retryable
                if attempt < MAX_ATTEMPTS - 1:
                    delay = RETRY_SLEEP_SECONDS[min(attempt, len(RETRY_SLEEP_SECONDS) - 1)]
                    print(
                        f"[!] {label} {batch_label}: attempt {attempt + 1}/{MAX_ATTEMPTS} "
                        f"failed ({exc}); retrying in {delay}s",
                        flush=True,
                    )
                    time.sleep(delay)
                else:
                    print(f"[!] {label} {batch_label}: giving up ({exc})", flush=True)

        if len(chunk) > MIN_SPLIT_SIZE:
            middle = len(chunk) // 2
            print(f"[!] {label} {batch_label}: splitting {len(chunk)} ids into "
                  f"{middle}+{len(chunk) - middle} and retrying", flush=True)
            process(chunk[:middle], batch_label + "a")
            process(chunk[middle:], batch_label + "b")
        else:
            failed.extend(chunk)
            print(f"[!] {label} {batch_label}: {len(chunk)} ids left unresolved",
                  flush=True)

    for index, start in enumerate(range(0, total, batch_size), start=1):
        chunk = ids[start:start + batch_size]
        print(
            f"[{index}/{n_batches}] {label}: {len(chunk)} ids "
            f"({start + 1}-{start + len(chunk)} of {total})",
            flush=True,
        )
        process(chunk, f"batch {index}")

    return parsed, failed


def taxid_sort_key(taxid):
    """Sort numeric taxids numerically, anything else lexicographically."""
    try:
        return (0, int(taxid))
    except (TypeError, ValueError):
        return (1, str(taxid))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_all_accessions_list", required=True,
                        help="Path to input assemblies file (one accession per line)")
    parser.add_argument("--output_taxonomy_tsv", required=True,
                        help="Path to output taxonomy TSV file")
    args = parser.parse_args()

    if shutil.which("efetch") is None:
        sys.exit("ERROR: `efetch` (NCBI Entrez Direct) is not available in PATH.")

    print(f"START GETTING TAXONOMY INFO using {args.input_all_accessions_list}", flush=True)
    print(
        "Using your NCBI API KEY" if os.environ.get("NCBI_API_KEY")
        else "No NCBI_API_KEY found: using the public rate limits",
        flush=True,
    )

    # Read, deduplicate and sort the accessions.
    accessions = []
    seen = set()
    with open(args.input_all_accessions_list) as accessions_list:
        for line in accessions_list:
            accession = line.strip()
            if not accession or accession in seen:
                continue
            seen.add(accession)
            accessions.append(accession)
    accessions.sort()

    if not accessions:
        print("[!] The input list is empty; writing an empty taxonomy.tsv", flush=True)
        with open(args.output_taxonomy_tsv, "w") as tax_file:
            tax_file.write(HEADER)
        return

    print(f"[i] {len(accessions)} unique accessions to resolve", flush=True)

    # 1) assembly accession -> taxid (batched docsum).
    assembly_to_taxid, failed_assemblies = fetch_and_parse(
        accessions, "assembly", "docsum", parse_assembly_docsum, "assembly->taxid"
    )
    accession_to_taxid = {}
    for accession in accessions:
        taxid = (assembly_to_taxid.get(accession)
                 or assembly_to_taxid.get(base_accession(accession)))
        if taxid:
            accession_to_taxid[accession] = taxid
    missing_assemblies = [acc for acc in accessions if acc not in accession_to_taxid]
    print(f"[+] assembly -> taxid: {len(accession_to_taxid)}/{len(accessions)} resolved",
          flush=True)

    # 2) taxid -> lineage (batched taxonomy XML). Many assemblies share a taxid,
    # so only the unique taxids are queried.
    unique_taxids = sorted(
        {taxid for taxid in accession_to_taxid.values() if taxid and taxid != "NA"},
        key=taxid_sort_key,
    )
    print(f"[i] {len(unique_taxids)} unique taxids to resolve", flush=True)
    taxid_to_lineage, failed_taxids = fetch_and_parse(
        unique_taxids, "taxonomy", "xml", parse_taxonomy_xml, "taxid->lineage"
    )
    print(f"[+] taxid -> lineage: {len(taxid_to_lineage)}/{len(unique_taxids)} resolved",
          flush=True)

    # 3) Join and write the TSV (one header, deterministic order).
    unresolved = 0
    with open(args.output_taxonomy_tsv, "w") as tax_file:
        tax_file.write(HEADER)
        for accession in accessions:
            taxid = accession_to_taxid.get(accession)
            lineage = taxid_to_lineage.get(taxid) if taxid else None
            if not lineage:
                # No taxonomy for this assembly: keep an explicit NA row so the
                # accession is visible in the report instead of silently absent.
                tax_file.write(f"{accession}\tNA\tNA\tNA\n")
                unresolved += 1
                continue
            for lineage_taxid, scientific_name, rank in lineage:
                tax_file.write(
                    f"{accession}\t{lineage_taxid}\t{scientific_name}\t{rank}\n"
                )

    print(f"[=] Written {args.output_taxonomy_tsv}: {len(accessions)} assemblies, "
          f"{unresolved} without taxonomy", flush=True)
    if missing_assemblies:
        print(f"[!] {len(missing_assemblies)} accessions have no assembly record in "
              f"NCBI: {', '.join(missing_assemblies[:20])}"
              f"{' ...' if len(missing_assemblies) > 20 else ''}", flush=True)
    if failed_assemblies:
        print(f"[!] {len(failed_assemblies)} accessions failed during the lookup: "
              f"{', '.join(failed_assemblies[:20])}"
              f"{' ...' if len(failed_assemblies) > 20 else ''}", flush=True)
    if failed_taxids:
        print(f"[!] {len(failed_taxids)} taxids could not be resolved: "
              f"{', '.join(failed_taxids[:20])}"
              f"{' ...' if len(failed_taxids) > 20 else ''}", flush=True)


if __name__ == "__main__":
    main()
