#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"
WINDOW=2000

for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    echo "=== $SP ==="
    # Dosyayı bul, bulamazsa donmak yerine hata verip atla
    GTF_GZ=$(ls $BASE/$SP/gtf/*.gtf.gz 2>/dev/null)
    if [ -z "$GTF_GZ" ]; then
        echo "  ERROR: GTF file not found for $SP!
        continue
    fi
    
    OUT=$BASE/$SP/promoters
    mkdir -p $OUT

    if [[ "$SP" == "hg38" || "$SP" == "mm10" ]]; then
        # GENCODE format - ZCAT ile direkt RAM üzerinden akış (Stream)
        zcat $GTF_GZ | awk -v W=$WINDOW 'BEGIN{OFS="\t"}
        $3=="gene" && /protein_coding/ {
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="") next
            if ($7=="+") { tss=$4; start=tss-W; end=tss+W }
            else         { tss=$5; start=tss-W; end=tss+W }
            if (start<1) start=1
            print $1, start, end, gname, ".", $7
        }' | sort -k1,1 -k2,2n > $OUT/promoters_all.bed
    else
        # Ensembl format - ZCAT ile direkt RAM üzerinden akış
        zcat $GTF_GZ | awk -v W=$WINDOW 'BEGIN{OFS="\t"}
        $3=="gene" && /protein_coding/ {
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="") {
                match($0, /gene_id "([^"]+)"/, b); gname=b[1];
            }
            if (gname=="") next
            if ($7=="+") { tss=$4; start=tss-W; end=tss+W }
            else         { tss=$5; start=tss-W; end=tss+W }
            if (start<1) start=1
            print "chr"$1, start, end, gname, ".", $7
        }' | grep -v "^chrMT\|^chrKB\|^chrKZ\|^chrML\|^chrJH\|^chrGL" \
           | sort -k1,1 -k2,2n > $OUT/promoters_all.bed
    fi

    echo "  Promoters: $(wc -l < $OUT/promoters_all.bed)"
    echo "  Sample: $(head -2 $OUT/promoters_all.bed | tail -1)"
done

echo ""
echo "=== COMPLETED ==="
