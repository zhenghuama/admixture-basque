PREFIX=step1_admixture
SCRIPT_DIR="."
INTERMEDIATE_DIR="${SCRIPT_DIR}/intermediate"
mkdir -p "$INTERMEDIATE_DIR"

REF_RESULTS_DIR="${SCRIPT_DIR}/reference_results"
mkdir -p "$REF_RESULTS_DIR"

# Keep all preprocessing artifacts under step1/intermediate/
BASE_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}"
PRUNED_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}.pruned"
RESULTS_PREFIX="${REF_RESULTS_DIR}/${PREFIX}"

for K in 3 5
do
    #echo "Not implemented - ADMIXTURE for K=$K"
    admixture "${PRUNED_PREFIX}.bed" "$K" > "${RESULTS_PREFIX}.${K}.0"
done