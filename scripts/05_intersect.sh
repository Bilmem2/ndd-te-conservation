#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"
RESULTS="$REPO_ROOT/results"
mkdir -p $RESULTS

for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    echo "=== $SP ==="
    OUT=$RESULTS/$SP
    mkdir -p $OUT

    for CAT in HighConfNDD Housekeeping; do
        PROM=$BASE/$SP/promoters/promoters_${CAT}.bed

        # LINE-1: tüm türler
        bedtools intersect -a $PROM -b $BASE/$SP/rmsk/LINE1.bed -c \
            > $OUT/${CAT}_LINE1.bed
        echo "  $CAT LINE1: done"

        # Alu: sadece primatlar
        if [[ "$SP" != "mm10" && "$SP" != "canFam4" ]]; then
            bedtools intersect -a $PROM -b $BASE/$SP/rmsk/Alu.bed -c \
                > $OUT/${CAT}_Alu.bed
            echo "  $CAT Alu: done"
        fi

        # B1/B2: sadece fare (primat Alu'nun fonksiyonel SINE analogu)
        if [[ "$SP" == "mm10" ]]; then
            bedtools intersect -a $PROM -b $BASE/$SP/rmsk/B1B2.bed -c \
                > $OUT/${CAT}_B1B2.bed
            echo "  $CAT B1/B2: done"
        fi
    done
done

echo "=== COMPLETED ==="
