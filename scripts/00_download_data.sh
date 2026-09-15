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
# Each data directory is named for the UCSC assembly it holds. The mouse lemur
# is absent here: its RepeatMasker track comes from the GenArk hub instead.
echo "=== RepeatMasker (UCSC) ==="
for asm in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 saiBol1 mm10 canFam6; do
    echo "- $asm"
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/$asm/database/rmsk.txt.gz" \
          "$DATA/$asm/rmsk/rmsk.txt.gz"
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
  [calJac4]=callithrix_jacchus [canFam6]=canis_lupus_familiaris
)
for folder in ponAbe3 nomLeu3 rheMac10 calJac4 canFam6; do
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
echo "Genomes and annotations complete."

# ── Auxiliary resources (confounder controls, matched controls, figures) ──────
# Required by 12, 15, 16 and 27. URLs verified 2026-09-10.
echo ""
echo "=== Auxiliary resources ==="

# gnomAD v4.1 gene constraint (LOEUF = lof.oe_ci.upper, pLI = lof.pLI), ~95 MB.
# NB the bucket path is release/4.1/, not release/v4.1/.
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
fetch "$GNOMAD_URL" "$DATA/hg38/gnomad_constraint.tsv"
# 07_pli_correlation.py expects it one level up; link rather than fetch twice.
if [[ -s "$DATA/hg38/gnomad_constraint.tsv" && ! -s "$DATA/gnomad_constraint.tsv" ]]; then
    cp "$DATA/hg38/gnomad_constraint.tsv" "$DATA/gnomad_constraint.tsv"
    echo "  copy -> gnomad_constraint.tsv (path expected by 07)"
fi

# GTEx v8 median TPM by tissue, ~7 MB. Bucket moved to adult-gtex.
fetch "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz" \
      "$DATA/hg38/gtex_gene_median_tpm.gct.gz"

# ClinVar variant summary, ~442 MB. 10_cross_disease.py reads it from data/.
fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" \
      "$DATA/clinvar_variants.txt.gz"

# ENCODE SCREEN Registry V4 candidate cis-regulatory elements, ~129 MB.
fetch "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed" \
      "$DATA/hg38/GRCh38-cCREs.bed"

# UCSC CpG islands (cpgIslandExt). Downloaded gzipped, then cut to BED.
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz" \
      "$DATA/hg38/cpgIslandExt.txt.gz"
if [[ -s "$DATA/hg38/cpgIslandExt.txt.gz" && ! -s "$DATA/hg38/cpg_islands.bed" ]]; then
    echo "  build cpg_islands.bed"
    gunzip -c "$DATA/hg38/cpgIslandExt.txt.gz" \
        | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $5}' > "$DATA/hg38/cpg_islands.bed"
fi

# Beagle GRCh38 genetic maps (sex-averaged cM/Mb), ~47 MB zipped.
# 15_context_controls.py globs recomb/chr_in_chrom_field/plink.chrchr*.GRCh38.map
fetch "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" \
      "$DATA/hg38/recomb/plink.GRCh38.map.zip"
RECOMB_DIR="$DATA/hg38/recomb/chr_in_chrom_field"
if [[ -s "$DATA/hg38/recomb/plink.GRCh38.map.zip" ]] \
   && ! compgen -G "$RECOMB_DIR/plink.chrchr*.GRCh38.map" >/dev/null; then
    echo "  unpack genetic maps"
    mkdir -p "$RECOMB_DIR"
    if command -v unzip >/dev/null 2>&1; then
        unzip -o -q -j "$DATA/hg38/recomb/plink.GRCh38.map.zip" -d "$RECOMB_DIR"
    else
        python -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
               "$DATA/hg38/recomb/plink.GRCh38.map.zip" "$RECOMB_DIR"
    fi
    # the archive ships plink.chr<N>.GRCh38.map; the loader expects chrchr<N>
    for f in "$RECOMB_DIR"/plink.chr[0-9XY]*.GRCh38.map; do
        [[ -e "$f" ]] || continue
        base="$(basename "$f")"
        case "$base" in
            plink.chrchr*) continue ;;
        esac
        mv "$f" "$RECOMB_DIR/plink.chrchr${base#plink.chr}"
    done
fi

echo ""
echo "=== DONE (auxiliary) ==="

# ── STILL MANUAL ─────────────────────────────────────────────────────────────
# ENCODE CTCF / fetal-brain DNase peak sets used by 08_encode_overlap_v2.py are
# accession-specific and must be pulled from the ENCODE portal by hand:
#   ENCFF955AQD ENCFF631TDE ENCFF670PXX ENCFF667IEN ENCFF362PZG ENCFF016LYI
#   -> data/hg38/{ctcf_peaks.bed, brain_dnase.bed.gz}
# 08 is an exploratory script and is not required to reproduce the manuscript.

echo ""
echo "Next step: bash scripts/02_rmsk_to_bed.sh"
