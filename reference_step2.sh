PREFIX=step2_admixture
SCRIPT_DIR="."
INTERMEDIATE_DIR="${SCRIPT_DIR}/data_and_results/step2/intermediate"
mkdir -p "$INTERMEDIATE_DIR"

REF_RESULTS_DIR="${SCRIPT_DIR}/data_and_results/step2/reference_results"
mkdir -p "$REF_RESULTS_DIR"

# Keep all preprocessing artifacts under step2/intermediate/
BASE_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}"
PRUNED_PREFIX="${INTERMEDIATE_DIR}/${PREFIX}.pruned"
RESULTS_PREFIX="${REF_RESULTS_DIR}/${PREFIX}"

for K in 3 5
do
    #echo "Not implemented - ADMIXTURE for K=$K"
    admixture "${PRUNED_PREFIX}.bed" "$K" > "${RESULTS_PREFIX}.${K}.0"
done

for K in 3 5
do
    mv "${SCRIPT_DIR}/${PREFIX}.pruned.${K}.P" "${RESULTS_PREFIX}.${K}.P"
    mv "${SCRIPT_DIR}/${PREFIX}.pruned.${K}.Q" "${RESULTS_PREFIX}.${K}.Q"
done