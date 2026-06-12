#!/bin/bash
# 00_download_data.sh — fetch genome GTF + RepeatMasker tracks for all 7 assemblies
# into the directory layout expected by scripts/02–05.
#
# The multi-GB raw files are intentionally NOT stored in the repo; this script
# reconstructs data/<assembly>/{gtf,rmsk}/ from public sources:
#   • GTF        — GENCODE (human, mouse) and Ensembl release 112 (other species)
#   • RepeatMasker — UCSC goldenPath rmsk.txt.gz
# Decompression is handled downstream by 02_rmsk_to_bed.sh / 03_get_promoters.sh.
#
# Idempotent: existing non-empty files are skipped, so re-running resumes a
# partial download. See README "Data Sources" for the canonical source table.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA="$REPO_ROOT/data"

# fetch <url> <outfile> — skip if present & non-empty; atomic; retry; warn on fail.
fetch() {
    local url="$1" out="$2"
    if [[ -s "$out" ]]; then
        echo "  skip (exists): ${out#"$DATA"/}"
        return 0
    fi
    mkdir -p "$(dirname "$out")"
    echo "  GET ${out#"$DATA"/}"
    if command -v curl >/dev/null 2>&1; then
        curl -fSL --retry 3 --retry-delay 5 -o "$out.part" "$url" && mv "$out.part" "$out" \
            || { echo "  WARN: download failed: $url" >&2; rm -f "$out.part"; return 1; }
    else
        wget -O "$out.part" "$url" && mv "$out.part" "$out" \
            || { echo "  WARN: download failed: $url" >&2; rm -f "$out.part"; return 1; }
    fi
}

# ── RepeatMasker (UCSC goldenPath) ────────────────────────────────────────────
# folder -> UCSC assembly. NB: dog data lives under data/canFam4/ for historical
# reasons but is actually ROS_Cfam_1.0 == UCSC canFam6 (coordinate-compatible;
# see README note). The folder name is kept as-is to preserve all path references.
echo "=== RepeatMasker (UCSC) ==="
declare -A RMSK_ASM=(
  [hg38]=hg38 [ponAbe3]=ponAbe3 [nomLeu3]=nomLeu3 [rheMac10]=rheMac10
  [calJac4]=calJac4 [mm10]=mm10 [canFam4]=canFam6
)
for folder in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    asm="${RMSK_ASM[$folder]}"
    echo "- $folder (UCSC $asm)"
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/$asm/database/rmsk.txt.gz" \
          "$DATA/$folder/rmsk/rmsk.txt.gz"
done

# ── GTF: GENCODE (human, mouse) ───────────────────────────────────────────────
# Human filename MUST stay gencode.v47.gtf.gz — 09/10_*.py reference it directly.
echo "=== GTF: GENCODE ==="
echo "- hg38 (GENCODE v47)"
fetch "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz" \
      "$DATA/hg38/gtf/gencode.v47.gtf.gz"
echo "- mm10 (GENCODE vM25)"
fetch "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" \
      "$DATA/mm10/gtf/gencode.vM25.gtf.gz"

# ── GTF: Ensembl release 112 (other primates + dog) ───────────────────────────
# 03_get_promoters.sh globs data/<folder>/gtf/*.gtf.gz, so the exact filename is
# irrelevant; we auto-discover the primary-assembly *.112.chr.gtf.gz from the
# Ensembl FTP listing (falls back to the toplevel *.112.gtf.gz if no .chr build).
echo "=== GTF: Ensembl release 112 ==="
declare -A ENS_DIR=(
  [ponAbe3]=pongo_abelii [nomLeu3]=nomascus_leucogenys [rheMac10]=macaca_mulatta
  [calJac4]=callithrix_jacchus [canFam4]=canis_lupus_familiaris
)
for folder in ponAbe3 nomLeu3 rheMac10 calJac4 canFam4; do
    sp="${ENS_DIR[$folder]}"
    base="https://ftp.ensembl.org/pub/release-112/gtf/$sp/"
    echo "- $folder (Ensembl $sp)"
    listing=$(curl -fsSL "$base" 2>/dev/null)
    fname=$(echo "$listing" | grep -oE '[A-Za-z_]+\.[A-Za-z0-9_.-]+\.112\.chr\.gtf\.gz' | sort -u | head -1)
    [[ -z "$fname" ]] && fname=$(echo "$listing" | grep -oE '[A-Za-z_]+\.[A-Za-z0-9_.-]+\.112\.gtf\.gz' | sort -u | head -1)
    if [[ -z "$fname" ]]; then
        echo "  WARN: could not auto-discover GTF at $base — fetch manually." >&2
        continue
    fi
    fetch "$base$fname" "$DATA/$folder/gtf/$fname"
done

echo ""
echo "=== DONE ==="
echo "Raw data written under: $DATA/<assembly>/{gtf,rmsk}/"
echo "Next step: bash scripts/02_rmsk_to_bed.sh"

# ── OPTIONAL auxiliary data (exploratory / confounder scripts) ────────────────
# Not downloaded automatically (version-specific and/or require local processing).
# See README for sources:
#   • ClinVar  variant_summary.txt.gz  -> data/clinvar_variants.txt.gz   (10_cross_disease.py)
#   • UCSC CpG islands (cpgIslandExt)   -> data/hg38/cpg_islands.bed      (12_figures_final.py, Fig 5)
#   • gnomAD v4.1 constraint            -> data/gnomad_constraint.tsv     (07_pli_correlation.py)
#   • ENCODE CTCF / brain DNase         -> data/hg38/{ctcf_peaks.bed,brain_dnase.bed.gz} (08_encode_overlap_v2.py)
