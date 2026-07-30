#!/bin/bash
# 00_download_data.sh — fetch genome GTF + RepeatMasker tracks for all 9 assemblies
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
  [calJac4]=calJac4 [saiBol1]=saiBol1 [mm10]=mm10 [canFam4]=canFam6
)
for folder in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 saiBol1 mm10 canFam4; do
    asm="${RMSK_ASM[$folder]}"
    echo "- $folder (UCSC $asm)"
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/$asm/database/rmsk.txt.gz" \
          "$DATA/$folder/rmsk/rmsk.txt.gz"
done

# ── hg38 genome sequence (UCSC bigZips) ───────────────────────────────────────
# Promoter GC content is read straight off the 2bit with twobitreader rather than
# via a FASTA, so this is a hard dependency of 18_matched_control.py,
# 19_rebaseline.py, 27_constraint_matched.py, 32_matching_sensitivity.py,
# 33_insertion_opportunity.py and 40_gc_analysis.py.
echo "=== hg38 genome sequence ==="
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit" \
      "$DATA/hg38/hg38.2bit"

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

# ── saiBol1 GTF (fixed output name) ───────────────────────────────────────────
# Handled outside the loop above: 39_squirrel_line1.py opens data/saiBol1/gtf/
# saiBol1.gtf.gz by name, and SaiBol1.0 has no chromosome-level (*.112.chr.gtf.gz)
# build, so the auto-discovery fallback would store a different filename.
echo "=== GTF: Ensembl release 112 (squirrel monkey) ==="
echo "- saiBol1 (Ensembl saimiri_boliviensis_boliviensis)"
fetch "https://ftp.ensembl.org/pub/release-112/gtf/saimiri_boliviensis_boliviensis/Saimiri_boliviensis_boliviensis.SaiBol1.0.112.gtf.gz" \
      "$DATA/saiBol1/gtf/saiBol1.gtf.gz"

# ── mmur3: mouse lemur (UCSC GenArk, flat layout) ─────────────────────────────
# Microcebus murinus has no UCSC goldenPath assembly, so its RepeatMasker track
# comes from the GenArk hub for GCF_000165445.2 (Mmur_3.0) as repeatMasker.out.gz
# rather than database/rmsk.txt.gz. Those coordinates carry RefSeq sequence names,
# so chromAlias.txt is required to map them onto the Ensembl annotation. The
# chromosome-level GTF (*.112.chr.gtf.gz) is the one used; 37_lemur_line1.py reads
# all three from a flat data/mmur3/ directory, not the gtf/ + rmsk/ split.
echo "=== mmur3: Ensembl GTF + UCSC GenArk RepeatMasker ==="
GENARK="https://hgdownload.soe.ucsc.edu/hubs/GCF/000/165/445/GCF_000165445.2"
echo "- mmur3 (Ensembl microcebus_murinus)"
fetch "https://ftp.ensembl.org/pub/release-112/gtf/microcebus_murinus/Microcebus_murinus.Mmur_3.0.112.chr.gtf.gz" \
      "$DATA/mmur3/mmur3.gtf.gz"
echo "- mmur3 (GenArk GCF_000165445.2)"
fetch "$GENARK/GCF_000165445.2.repeatMasker.out.gz" "$DATA/mmur3/repeatMasker.out.gz"
fetch "$GENARK/GCF_000165445.2.chromAlias.txt"      "$DATA/mmur3/chromAlias.txt"

echo ""
echo "=== DONE ==="
echo "Raw data written under: $DATA/<assembly>/{gtf,rmsk}/ (mmur3 is flat: data/mmur3/)"
echo "Next step: bash scripts/02_rmsk_to_bed.sh"

# ── OPTIONAL auxiliary data (exploratory / confounder scripts) ────────────────
# Not downloaded automatically (version-specific and/or require local processing).
# See README for sources:
#   • ClinVar  variant_summary.txt.gz  -> data/clinvar_variants.txt.gz   (10_cross_disease.py)
#   • UCSC CpG islands (cpgIslandExt)   -> data/hg38/cpg_islands.bed      (12_figures_final.py, Fig 5)
#   • gnomAD v4.1 constraint            -> data/gnomad_constraint.tsv     (07_pli_correlation.py)
#   • ENCODE CTCF / brain DNase         -> data/hg38/{ctcf_peaks.bed,brain_dnase.bed.gz} (08_encode_overlap_v2.py)
