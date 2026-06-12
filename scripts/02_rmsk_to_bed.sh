#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"

for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    echo "=== $SP ==="
    RMSK=$BASE/$SP/rmsk/rmsk.txt.gz
    OUT=$BASE/$SP/rmsk
    
    # Aç
    gunzip -k $RMSK 2>/dev/null || true
    RMSK_TXT=$BASE/$SP/rmsk/rmsk.txt

    # LINE-1 (tüm türler)
    awk 'NR>1 && $12=="LINE" && $13=="L1" {
        printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
    }' $RMSK_TXT | sort -k1,1 -k2,2n > $OUT/LINE1.bed
    echo "  LINE1: $(wc -l < $OUT/LINE1.bed) element"

    # Alu (primatlar: hg38, rheMac10, calJac4, ponAbe3)
    if [[ "$SP" == "hg38" || "$SP" == "ponAbe3" || "$SP" == "nomLeu3" || "$SP" == "rheMac10" || "$SP" == "calJac4" ]]; then
        awk 'NR>1 && $12=="SINE" && $13=="Alu" {
            printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
        }' $RMSK_TXT | sort -k1,1 -k2,2n > $OUT/Alu.bed
        echo "  Alu: $(wc -l < $OUT/Alu.bed) element"
    fi

done

echo ""
echo "=== TAMAMLANDI ==="
