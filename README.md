# Inference of within-host viral dynamics from cross-sectional viral load distributions

Code for Green, Cori and Ferguson, *Inference of within-host viral dynamics from
cross-sectional population-level viral load distributions: application to SARS-CoV-2 in the UK*.

The model infers viral growth rate, decay rate, peak viral load and the delay to symptom onset,
each with its between-individual distribution, from repeated cross-sectional viral load
measurements rather than from followed-up individuals. It is applied to Pillar 2 community
testing data from England for the wild-type, Alpha, Delta and Omicron variants.

## Requirements

R 4.1+, and a C++11 toolchain with OpenMP (Rtools on Windows; `libomp` on macOS, since Apple
Clang does not ship it).

```r
install.packages(c("Rcpp", "dqrng", "BH", "sitmo",
                   "TruncatedNormal", "mvtnorm", "truncnorm"))
install.packages(c("rstan", "dplyr", "ggplot2", "matrixStats"))  # longitudinal fits
install.packages(c("posterior", "ggh4x", "patchwork"))           # diagnostics and figures
```

## Running it

```
Rscript run_all.R check          environment, kernels and data        ~2 min
Rscript run_all.R demo           one short fit, end to end            ~5 min
Rscript run_all.R simulations    the 24 simulation-recovery fits      ~1 core-day each
Rscript run_all.R score          recovery against the known truth
Rscript run_all.R figures        the simulation-study figures
Rscript run_all.R longitudinal   refit the viral kinetics priors      a few hours
Rscript run_all.R regenerate     rebuild the synthetic datasets
```

Start with `check`. It compiles the kernels and confirms the datasets are intact.

Each simulation fit is 250,000 iterations and about a core-day, so running all 24 sequentially is
impractical. Launch them as separate processes instead, one per core:

```bash
FIT_FORM=incidence SYMP_DELAY_LIM=28 FIT_N_ITER=250000 \
FIT_COV_START=80000,100000,120000,140000 FIT_COV_WINDOW=20000 \
Rscript scripts/fit_model.R A0 1 1
```

Chains checkpoint to `outputs/` every 5,000 iterations. Posterior summaries use the draws after
the final covariance adaptation, recorded on each fit as `x$analysis_from`.

## What needs the Pillar 2 data

The variant fits and eleven of the eighteen figures need the Pillar 2 testing data, which is
UKHSA data held under agreement and is not ours to release. Access can be requested from UKHSA
through the [Office for Data Release](https://www.gov.uk/government/publications/accessing-ukhsa-protected-data)
or at `DataAccess@ukhsa.gov.uk`.

The code is the same either way. Fitting a variant is

```bash
FIT_FORM=incidence SYMP_DELAY_LIM=28 Rscript scripts/fit_model.R wt 1 1
```

and differs from a synthetic scenario only in which array it reads. With the arrays in `data/`,
`figure3()`, `figure4()`, `figure5()`, `figure4_priors()`, `figure_corr()`, `figure_pairs()`,
`figure_prior_posterior()` and `figure_phase_ppc()` in `scripts/paper_outputs.R` run unchanged,
as do Tables 1, S1, S2 and S6.

Everything else runs here: the simulation study, its four figures, and the longitudinal priors.

## Synthetic datasets

Twenty-four datasets in `data/`, covering two viral kinetics profiles by three positions of
symptom onset relative to peak by four sample sizes.

```
A  less variable, onset mostly before peak      0  ~1,000 positive tests
B  less variable, onset about half before       1  ~10,000
C  less variable, onset mostly after peak       2  ~100,000
D  more variable, onset mostly before peak      3  ~1,000,000
E  more variable, onset about half before
F  more variable, onset mostly after peak
```

So `C2` is the less variable profile with late onset at about 100,000 tests. Each file is
`list(out, pars_test)`: a 22 x 105 x 7 array of counts by viral load bin, calendar day and days
from onset to test, plus the parameter vector that generated it. `input_generator` surfaces that
vector as `x$truth`, which is how recovery is scored.

The M-R arrays hold the same generating vectors under a 20-day threshold-to-onset window, where
A-F use the 28-day window the paper reports. `run_all.R regenerate` rebuilds A-F from them and
checks the result byte for byte.

Kinetics priors are widened twofold for synthetic targets and the spline coefficients given a
coefficient of variation of one, so recovery tests what the data carries rather than what the
prior already knew. Generation uses `test_pop = 1e7` against the fits' `1e6`, so the Monte Carlo
error in the data is independent of that in the likelihood.

## Longitudinal priors

The viral kinetics priors come from a hierarchical model fitted to two cohorts, both published by
the groups that produced them:

- `data/ct_data_clean.csv` — NBA, Kissler et al.,
  [gradlab/CtTrajectories](https://github.com/gradlab/CtTrajectories). Viral loads are used as
  given; no Ct conversion is applied, as the published values are already calibrated.
- `data/20210920_ATACCC_L.csv` — ATACCC, Singanayagam et al.,
  [mrc-ide/ATACCC-kinetic-model-and-data](https://github.com/mrc-ide/ATACCC-kinetic-model-and-data).
  Ct is converted with the study's own relationship, `log10(133.33 * exp((37.933 - CtT1) / 1.418))`.

Please cite the original studies if you use either. Pillar 2 uses the García-Fiñana TaqPath
calibration, so three platform-appropriate conversions are used across the project.

`run_all.R longitudinal` writes `long_post_6.rds` and `long_cov_6.rds`, the prior mean and
covariance the main model reads, and the ATACCC-only equivalents used as a sensitivity check.
The `6` is a trajectory-form index passed as `formula = 6`, not a file version; there is one Stan
model. Twenty-nine participants are excluded through `reject_IDS`, being those whose series does
not identify a trajectory.

Both summaries are committed, so the simulation study runs without refitting.

## Layout

```
run_all.R                      every stage
scripts/
  fit_model.R                  fits one chain to one target
  model.R                      likelihood, priors, input assembly
  sampler.R                    Metropolis-Hastings
  simulate_incidence.cpp       trajectories, incidence indexing (reported model)
  simulate_onset.cpp           trajectories, onset indexing (earlier formulation)
  score_recovery.R             recovery against known parameters
  make_figures.R               the figures buildable here
  paper_outputs.R              all figures and tables
  make_scenarios_incidence.R   builds the synthetic datasets
  synthetic_truth.R            generating-parameter snapshot
  longitudinal_fit.R           the Stan prior fit
data/                          synthetic datasets and the two cohorts
longitudinal_fits/             Stan model and the prior summaries
outputs/                       chains land here, not committed
```

`run_all.R` calls `fit_model.R`, which sources `model.R` and `sampler.R`; `sampler.R` compiles
both kernels. To change the model, edit `model.R`; to change how it is sampled, `sampler.R`.
`FIT_FORM` selects which kernel is used.

## Reproducibility

The posterior defined here is the one the reported fits sampled: re-evaluating the log posterior
at states those chains visited reproduces their recorded values to within 1e-12.

Fresh runs will not return identical chains. The reported fits took their proposal covariance
from preliminary runs, so several are continuations rather than single invocations and their
lengths differ. The kernel also allocates random number streams per core, so core count changes
the Monte Carlo realisation of the likelihood. Expect the coverage count near the published 230
of 240 rather than exactly on it.

## Data availability

The Pillar 2 data is confidential and cannot be made publicly available; see above for the
request route. The synthetic datasets are provided so that the method and the simulation study
can be reproduced without it.

## Licence

MIT, see `LICENSE`.
