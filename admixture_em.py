"""
Method 1:
Naive E-M implementation of Admixture, no further optimizations.
Slow convergence.
Achieve a reasonable result, final log likelihood is worse than Admixture

Pre-requisites:
Use Plink to perform pre-processing, inclusing LD pruning

Reads PLINK binary files from preprocessing output:
  - intermediate/step1_admixture.pruned.{bed,bim,fam}) 
and writes:
  - results/step1_admixture.pruned.{K}.Q
  - results/step1_admixture.pruned.{K}.P
  - results/step1_admixture.{K}.0 (log file)

"""

import time
from pathlib import Path
from typing import Callable
import numpy as np


def count_lines(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def read_plink_bed(prefix: Path) -> np.ndarray:
    """
    Read SNP-major PLINK genotype binary.

    Returns:
        Genotype matrix (N samples x M SNPs)
        int8 matrix with genotype counts for allele1 in .bim:
          0: homozygous allele1
          1: heterozygous
          2: homozygous allele2
         -1: missing / error
    """
    fam_path = Path(str(prefix) + ".fam") # Sample list, 1 line / ppl
    bim_path = Path(str(prefix) + ".bim") # SNP list, 1 line / variant
    bed_path = Path(str(prefix) + ".bed") # binary, SNP major, 2 bits / genotype

    n_samples = count_lines(fam_path)
    n_snps = count_lines(bim_path)
    bytes_per_snp = (n_samples + 3) // 4

    raw = np.fromfile(bed_path, dtype=np.uint8)
    # PLINK .bed files start with 3-byte header: 0x6C 0x1B 0x01 for SNP-major.
    expected_size = 3 + n_snps * bytes_per_snp
    if raw.size != expected_size:
        raise ValueError(
            f"{bed_path} size mismatch: expected {expected_size} bytes, found {raw.size}."
        )

    blocks = raw[3:].reshape(n_snps, bytes_per_snp)

    # Unpack 2-bit genotype codes (low bits first, 4 samples per byte).
    codes = np.empty((n_snps, bytes_per_snp * 4), dtype=np.uint8)
    codes[:, 0::4] = blocks & 0b11
    codes[:, 1::4] = (blocks >> 2) & 0b11
    codes[:, 2::4] = (blocks >> 4) & 0b11
    codes[:, 3::4] = (blocks >> 6) & 0b11
    codes = codes[:, :n_samples].T  # N x M

    # PLINK mapping:
    # 00: hom allele1, 10: het, 11: hom allele2, 01: missing.
    genotype = np.full((n_samples, n_snps), -1, dtype=np.int8)
    genotype[codes == 0b00] = 0
    genotype[codes == 0b10] = 1
    genotype[codes == 0b11] = 2
    return genotype


def initialize_qf(
    genotype: np.ndarray,
    geno_valid_mask: np.ndarray,
    k: int,
    rng: np.random.Generator,
    min_prob: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Initialize Q and F matrices for EM.
    Inputs:
    - genotype: N x M genotype counts in {0,1,2}, with -1 for missing (floating point at this step).
    - geno_valid_mask: N x M boolean mask indicating valid genotypes.
    - k: Number of ancestral populations.
    - rng: Random number generator.
    - min_prob: Minimum probability for numerical stability.
    Outputs:
    - q_init: N x K matrix of initial admixture proportions.
    - f_init: K x M matrix of initial allele frequencies.
    """
    genotype = np.where(geno_valid_mask > 0, genotype, 0.0)

    n_samples, n_snps = genotype.shape

    # Initialize Q matrix (admixture proportions) randomly.
    # Sum k q_ik should be 1 for each role (individual) 
    q_init = rng.random((n_samples, k))
    q_init /= q_init.sum(axis=1, keepdims=True)

    # Initialize F matrix empirical allele frequency, with small noise per cluster.
    called = geno_valid_mask.sum(axis=0).clip(min=1)
    print("Genotype summary:")
    print("Samples:", n_samples, "SNPs:", n_snps)
    print(called.min(), called.mean(), called.max())
    allel_prob = (genotype * geno_valid_mask).sum(axis=0) / (2.0 * called)
    f_init = np.repeat(allel_prob.reshape(1, -1), k, axis=0)
    f_init += rng.normal(0.0, 0.01, size=f_init.shape)
    print("Initial F summary:")
    print(allel_prob.min(), allel_prob.mean(), allel_prob.max())
    print(f_init.min(), f_init.mean(), f_init.max())
    np.clip(f_init, min_prob, 1.0 - min_prob, out=f_init)
    return q_init, f_init


def em_admixture(
    genotype: np.ndarray,
    k: int,
    max_iters: int,
    early_stop_rate: float,
    random_seed: int,
    min_prob: float,
    iteration_logger: Callable[[int, float, float], None] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Baseline EM for admixture model.
    Inputs:
        - genotype: N x M genotype counts in {0,1,2}, with -1 for missing.
        - k: Number of ancestral populations.
        - max_iters: Maximum number of EM iterations.
        - early_stop_rate: Relative log-likelihood improvement threshold for early stopping.
        - random_seed: 42 lol
        - min_prob: Minimum probability for numerical stability.
        - iteration_logger: Callback for logging iteration info (iteration, log-likelihood, delta
    Outputs:
        - q: N x K matrix of admixture proportions.
        - f: K x M matrix of allele frequencies.
    """
    genotype = genotype.astype(np.float64, copy=False)
    geno_valid_mask = (genotype >= 0).astype(np.float64)

    rng = np.random.default_rng(random_seed)

    # Q shape: N x K, F shape: K x M
    q, f = initialize_qf(
        genotype=genotype, 
        geno_valid_mask=geno_valid_mask, 
        k=k, 
        rng=rng, 
        min_prob=min_prob
    )

    llh_log = np.zeros(max_iters, dtype=np.float64)

    #logs: list[tuple[int, float, float]] = []
    prev_ll = -np.inf

    for it in range(1, max_iters + 1):
        p = q @ f
        np.clip(p, min_prob, 1.0 - min_prob, out=p)

        ll = np.sum(
            geno_valid_mask * (genotype * np.log(p) + (2.0 - genotype) * np.log1p(-p)) 
            # using log1p for numerical stability when p is close to 1.0
        )
        delta = ll - prev_ll if np.isfinite(prev_ll) else np.nan 
        llh_log[it - 1] = ll.astype(np.float64)

        if iteration_logger is not None:
            iteration_logger(it, ll, delta)

        if it > 1 and abs(delta) <= early_stop_rate * (1.0 + abs(prev_ll)):
            llh_log = llh_log[:it]
            break

        q_num = np.zeros_like(q)
        f_num = np.zeros_like(f)
        f_den = np.zeros_like(f)

        # E-step expected allele ancestry counts + M-step closed forms.
        for kk in range(k):
            qk = q[:, [kk]]          # N x 1
            fk = f[[kk], :]          # 1 x M

            alt_prob_kk = qk @ fk / p
            ref_prob_kk = qk @ (1.0 - fk) / (1.0 - p)
            
            alt_exp_kk = geno_valid_mask * genotype * alt_prob_kk
            ref_exp_kk = geno_valid_mask * (2.0 - genotype) * ref_prob_kk
            tot_exp_kk = alt_exp_kk + ref_exp_kk

            q_num[:, kk] = tot_exp_kk.sum(axis=1)
            f_num[kk, :] = alt_exp_kk.sum(axis=0)
            f_den[kk, :] = tot_exp_kk.sum(axis=0)

        # Update Q
        # Theoreetically Q_row_sum should be 2*n_snps, but we use the sum of expected q_num to be more robust to missing data
        q_row_sums = q_num.sum(axis=1, keepdims=True)
        q = q_num / q_row_sums

        # Update F.
        np.divide(f_num, f_den, out=f, where=f_den > 0.0)
        np.clip(f, min_prob, 1.0 - min_prob, out=f)

        prev_ll = ll
          # trim log to actual iterations

    return q, f, llh_log


def write_outputs(
    q: np.ndarray,
    f: np.ndarray,
    llh_log: np.ndarray,
    out_dir: Path,
    pruned_name: str,
    k: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    q_path = out_dir / f"{pruned_name}.{k}.Q"
    p_path = out_dir / f"{pruned_name}.{k}.P"
    llh_path = out_dir / f"{pruned_name}.{k}.npy"

    np.savetxt(q_path, q, fmt="%.6f")
    # .Q saves F matrix in SNP major M x K
    np.savetxt(p_path, f.T, fmt="%.6f")
    np.save(llh_path, llh_log)
    return None


def run_admixture(
    input_prefix: str,
    output_folder: str = "results",
    k_values: list[int] = [3, 5],
    max_iters: int = 1000,
    early_stop_rate: float = 1e-4,
    random_seed: int = 42,
    min_prob: float = 1e-5
) -> None:
    """
    Run baseline admixture EM for one or more K values.

    Designed for direct notebook usage:
      results = run_admixture(...)
    """
    input_prefix = Path(input_prefix)
    script_dir = Path(__file__).resolve().parent
    out_dir = Path(output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    genotype = read_plink_bed(input_prefix)

    pruned_name = input_prefix.name  # e.g. step1_admixture.pruned
    log_name = pruned_name[:-7] 

    for idx, k in enumerate(k_values):
        log_path = out_dir / f"{log_name}.{k}.0"
        start = time.time()
        with log_path.open("w") as log_handle:
            log_handle.write("****      Baseline ADMIXTURE: E-M      ****\n")
            log_handle.write(f"K={k}\n")
            log_handle.write(f"Samples={genotype.shape[0]} SNPs={genotype.shape[1]}\n\n")

            def _iter_logger(it: int, ll: float, delta: float) -> None:
                if np.isnan(delta):
                    log_handle.write(f"{it}\tLoglikelihood: {ll:.6f}\t(delta): NA\n")
                else:
                    log_handle.write(
                        f"{it}\tLoglikelihood: {ll:.6f}\t(delta): {delta:.6f}\n"
                    )
                log_handle.flush()

            q, f, llh_log = em_admixture(
                genotype=genotype,
                k=k,
                max_iters=max_iters,
                early_stop_rate=early_stop_rate,
                random_seed=random_seed + idx,
                min_prob=min_prob,
                iteration_logger=_iter_logger,
            )
            runtime_s = time.time() - start
            log_handle.write(f"\nRuntime(s): {runtime_s:.2f}\n")
            log_handle.flush()

        write_outputs(
            q=q,
            f=f,
            llh_log=llh_log,
            out_dir=out_dir,
            pruned_name=pruned_name,
            k=k,
        )


    return None
