#!/bin/bash
# NOTE: Archived/exploratory — superseded by scripts/03_get_promoters.sh.
# Retained for provenance only; not part of the main pipeline.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"
WINDOW=2000

echo "Promotor koordinatları çıkarılıyor (TSS ± ${WINDOW} bp, protein_coding only)..."

for SP in hg38 rheMac10 calJac4 mm10; do
    echo ""
    echo "=== $SP ==="
    
    GTF_GZ=$(ls $BASE/$SP/gtf/*.gtf.gz)
    OUT=$BASE/$SP/promoters
    GTF=$(echo $GTF_GZ | sed 's/.gz//')
    
    if [ "$SP" == "hg38" ] || [ "$SP" == "mm10" ]; then
        # GENCODE: gene_name var, protein_coding filtrele
        awk -v W=$WINDOW 'BEGIN{OFS="\t"}
        $3=="gene" && /protein_coding/ {
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="") next
            if ($7=="+") { tss=$4; start=tss-W; end=tss+W }
            else          { tss=$5; start=tss-W; end=tss+W }
            if (start<1) start=1
            print $1, start, end, gname, ".", $7
        }' $GTF | sort -k1,1 -k2,2n > $OUT/promoters_all.bed

    else
        # Ensembl: gene_name yoksa gene_id kullan, chr prefix ekle
        awk -v W=$WINDOW 'BEGIN{OFS="\t"}
        $3=="gene" && /protein_coding/ {
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="") {
                match($0, /gene_id "([^"]+)"/, b); gname=b[1];
            }
            if (gname=="") next
            if ($7=="+") { tss=$4; start=tss-W; end=tss+W }
            else          { tss=$5; start=tss-W; end=tss+W }
            if (start<1) start=1
            print "chr"$1, start, end, gname, ".", $7
        }' $GTF | grep -v "^chrMT\|^chrKB\|^chrKZ\|^chrML" \
               | sort -k1,1 -k2,2n > $OUT/promoters_all.bed
    fi
    
    echo "  Protein-coding promotor: $(wc -l < $OUT/promoters_all.bed)"
    echo "  Örnek: $(head -2 $OUT/promoters_all.bed | tail -1)"
    echo "  Boş gen adı: $(awk '$4==""' $OUT/promoters_all.bed | wc -l)"
done

echo ""
echo "=== TAMAMLANDI ==="
