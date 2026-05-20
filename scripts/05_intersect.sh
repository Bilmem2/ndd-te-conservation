#!/bin/bash
BASE=~/comparative_TE_study/data
RESULTS=~/comparative_TE_study/results
mkdir -p $RESULTS

for SP in hg38 rheMac10 calJac4 mm10 ponAbe3 canFam4; do
    echo "=== $SP ==="
    OUT=$RESULTS/$SP
    mkdir -p $OUT

    for CAT in HighConfNDD BroadNDD Housekeeping; do
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
    done
done

echo "=== TAMAMLANDI ==="
