#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"

for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 mm10 canFam4; do
    echo "=== $SP ==="
    RMSK=$BASE/$SP/rmsk/rmsk.txt.gz
    OUT=$BASE/$SP/rmsk

    if [[ ! -s "$RMSK" ]]; then
        echo "  ERROR: $RMSK not found — run 00_download_data.sh first." >&2
        continue
    fi

    # rmsk.txt.gz akıtılarak okunur; açılmış kopya (tür başına ~1.5 GB) yazılmaz.
    # LINE-1 (tüm türler)
    zcat $RMSK | awk 'NR>1 && $12=="LINE" && $13=="L1" {
        printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
    }' | sort -k1,1 -k2,2n > $OUT/LINE1.bed
    echo "  LINE1: $(wc -l < $OUT/LINE1.bed) element"

    # Alu (primatlar: hg38, rheMac10, calJac4, ponAbe3)
    if [[ "$SP" == "hg38" || "$SP" == "ponAbe3" || "$SP" == "nomLeu3" || "$SP" == "rheMac10" || "$SP" == "calJac4" ]]; then
        zcat $RMSK | awk 'NR>1 && $12=="SINE" && $13=="Alu" {
            printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
        }' | sort -k1,1 -k2,2n > $OUT/Alu.bed
        echo "  Alu: $(wc -l < $OUT/Alu.bed) element"
    fi

    # B1/B2 (fare): primat Alu'nun fonksiyonel SINE analogları.
    # RepeatMasker'da B1 SINE/Alu, B2 ise SINE/B2 ailesinde sınıflanır.
    if [[ "$SP" == "mm10" ]]; then
        zcat $RMSK | awk 'NR>1 && $12=="SINE" && ($13=="Alu" || $13=="B2") {
            printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
        }' | sort -k1,1 -k2,2n > $OUT/B1B2.bed
        echo "  B1/B2: $(wc -l < $OUT/B1B2.bed) element"
    fi

done

echo ""
echo "=== TAMAMLANDI ==="
