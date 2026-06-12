#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"
LISTS=$BASE/gene_lists

for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    echo "=== $SP ==="
    OUT=$BASE/$SP/promoters
    ALL=$OUT/promoters_all.bed

    # Mouse gen sembolleri ilk harf büyük — uppercase yap
    if [ "$SP" == "mm10" ]; then
        awk 'BEGIN{OFS="\t"} {$4=toupper($4); print}' $ALL > $OUT/promoters_upper.bed
        MATCH=$OUT/promoters_upper.bed
    else
        MATCH=$ALL
    fi

    for CAT in HighConfNDD BroadNDD Housekeeping; do
        grep -Fw -f $LISTS/${CAT}_genes.txt $MATCH > $OUT/promoters_${CAT}.bed
        echo "  $CAT: $(wc -l < $OUT/promoters_${CAT}.bed) gen"
    done
done

echo ""
echo "=== TAMAMLANDI ==="
