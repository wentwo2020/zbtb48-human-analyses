#!/usr/bin/env bash
# Run every analysis that needs only the base requirements.
# Scripts 04 and 05 are excluded: 04 needs requirements-singlecell.txt and 05
# downloads about 0.65 GB. Run those two explicitly when you want them.
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_gtex_age_association \
         02_gtex_reference_genes_and_composition \
         03_gtex_qtl_lead_variant \
         06_finngen_mr_and_enrichment; do
  echo "=== $s"
  python "src/$s.py"
done
echo
echo "All analyses completed and all assertions passed. Results in results/."
