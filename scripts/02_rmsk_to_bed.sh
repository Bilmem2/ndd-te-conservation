#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE="$REPO_ROOT/data"

# saiBol1 ships the standard UCSC rmsk table and 00_download_data.sh fetches
# it, but it was missing from this loop, so its Alu/LINE-1 beds were never
# derived and 30_squirrel_monkey.py had to parse the raw table itself.
for SP in hg38 ponAbe3 nomLeu3 rheMac10 calJac4 saiBol1 mm10 canFam6; do
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
    if [[ "$SP" == "hg38" || "$SP" == "ponAbe3" || "$SP" == "nomLeu3" || "$SP" == "rheMac10" || "$SP" == "calJac4" || "$SP" == "saiBol1" ]]; then
        zcat $RMSK | awk 'NR>1 && $12=="SINE" && $13=="Alu" {
            printf "%s\t%d\t%d\t%s\n",$6,$7,$8,$11
        }' | sort -k1,1 -k2,2n > $OUT/Alu.bed
        echo "  Alu: $(wc -l < $OUT/Alu.bed) element"

        # Nine of the human analyses (19, 25, 26, 27, 28, 32, 33, 36 and
        # probe_dosage) read the same records under the older flat name
        # data/hg38/alu_rmsk.bed. Nothing produced it, so a clean clone could
        # not run any of them. They all take columns 1-3 only, so a copy of
        # Alu.bed serves; kept as a copy rather than a link for Windows.
        if [[ "$SP" == "hg38" ]]; then
            cp "$OUT/Alu.bed" "$BASE/$SP/alu_rmsk.bed"
            echo "  alu_rmsk.bed: flat copy for scripts 19-36"
        fi
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
