PREFIX=step3_admixture
SCRIPT_DIR="."
INTERMEDIATE_DIR="${SCRIPT_DIR}/data_and_results/step3/intermediate"
mkdir -p "$INTERMEDIATE_DIR"

REF_RESULTS_DIR="${SCRIPT_DIR}/data_and_results/step3/reference_results"
mkdir -p "$REF_RESULTS_DIR"

# Keep all preprocessing artifacts under step3/intermediate/
BASE_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}"
PRUNED_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}.pruned"
RESULTS_PREFIX="${REF_RESULTS_DIR}/${PREFIX}"

for K in 3
do
    #echo "Not implemented - ADMIXTURE for K=$K"
    admixture "${PRUNED_PREFIX}.bed" "$K" > "${RESULTS_PREFIX}.${K}.0"
done

for K in 3
do
    mv "${SCRIPT_DIR}/${PREFIX}.pruned.${K}.P" "${RESULTS_PREFIX}.pruned.${K}.P"
    mv "${SCRIPT_DIR}/${PREFIX}.pruned.${K}.Q" "${RESULTS_PREFIX}.pruned.${K}.Q"
done