#!/bin/bash
BASE=~/comparative_TE_study/data
WINDOW=2000

for SP in hg38 rheMac10 calJac4 mm10 ponAbe3 canFam4; do
    echo "=== $SP ==="
    GTF_GZ=$(ls $BASE/$SP/gtf/*.gtf.gz)
    GTF=$(echo $GTF_GZ | sed 's/.gz//')
    OUT=$BASE/$SP/promoters
    gunzip -k $GTF_GZ 2>/dev/null || true

    if [[ "$SP" == "hg38" || "$SP" == "mm10" ]]; then
        # GENCODE format
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
        # Ensembl format — chr prefix ekle
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
        }' $GTF | grep -v "^chrMT\|^chrKB\|^chrKZ\|^chrML\|^chrJH\|^chrGL" \
               | sort -k1,1 -k2,2n > $OUT/promoters_all.bed
    fi

    echo "  Promotor: $(wc -l < $OUT/promoters_all.bed)"
    echo "  Örnek: $(head -2 $OUT/promoters_all.bed | tail -1)"
done

echo ""
echo "=== TAMAMLANDI ==="
