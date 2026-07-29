# Cross-disease results

`cross_disease_recomputed.csv` holds the numbers reported in the manuscript.

The earlier outputs of `10_cross_disease.py` (`cross_disease_results.csv`,
`*_promoters.bed`) were removed because they were produced before a parsing bug
was fixed. ClinVar reports several genes per variant in one semicolon-delimited
field and writes a bare `-` where no gene is assigned; selecting promoters with
`grep -Fw` let that `-` match the strand column of every minus-strand gene, so
the "Mendelian" set was roughly three-quarters arbitrary genes (9,912 minus- vs
2,704 plus-strand promoters, where a real gene set is near 50/50).

Both `10_cross_disease.py` and `26_cross_disease_recompute.py` now split the
field and match on exact gene symbol. Run either to regenerate the corrected
sets; `26_` works from the committed gene lists and does not need the ClinVar
download.
