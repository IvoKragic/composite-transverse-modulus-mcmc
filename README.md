# composite-transverse-modulus-mcmc

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23774307.svg)](https://doi.org/10.5281/zenodo.23774307)

Archived version of this release: https://doi.org/10.5281/zenodo.23774307

MATLAB code — Bayesian model selection of the transverse modulus of a composite ply

Code accompanying the manuscript

> **Micromechanical model selection and parameter estimation of the transverse elastic modulus
> of a unidirectional composite ply using a Bayesian approach**
> Stipe Perišić, Ivo Marinić-Kragić\*, Frane Vlak, Jani Barle
> *Journal of Composite Materials* (submitted)

It contains the MATLAB implementation of the Bayesian analysis (Metropolis–Hastings MCMC and
model-switching RJMCMC) for three micromechanical models of the transverse elastic modulus
*E₂* of a unidirectional ply: Halpin–Tsai with ξ = 2, Halpin–Tsai with the Hewitt–de Malherbe
reinforcement factor, and Puck.

## Requirements

- MATLAB R2026b (tested) or newer
- Statistics and Machine Learning Toolbox (`normpdf`, `normrnd`, `quantile`, `nlinfit`, `statset`)
- Curve Fitting Toolbox — only for the `fit`-based comparison in `mcmc/`

No other dependencies.

## Quick start

```matlab
cd release            % folder containing startup.m
startup               % adds models/, mcmc/, rjmcmc/, validation/, data/ to the path

run_all               % reproduces Table 6 and Supplementary Tables S1-S4 (quick settings)
run_all('full')       % same with the settings used in the paper (slower, ~15-25 min)
```

Each block of `run_all` prints a table that can be compared directly with the manuscript.

## Contents

```
models/       E_2HalpinTsai.m     Halpin–Tsai model ('original', 'Hewitt', 'Giner')
              E_2Puck.m           Puck's model (νm = 0.35 by default)

mcmc/         likelihoodHT.m      Gaussian log-likelihood
              sampleHT.m          single-model Metropolis sampler (as used for Figures 3–5)
              stipePlotHT.m       trace/histogram plots

rjmcmc/       likelihoodHTrjmcmc.m   log-likelihood used by the model-jumping chain
              sampleHTrjmcmc.m      within-model Metropolis step
              composite_bayes_rjmcmc_ver3.m  original driver (model selection, Figures 6–9)

validation/   rw_mcmc.m           Metropolis sampler with positive-support prior on σ
              ess_acf.m           effective sample size (initial positive sequence)
              rjmcmc_quick.m      compact model-switching chain
              ess/replication/sensitivity scripts used for Tables S1–S4

data/         E2podaci.dat        15 experimental (Vf, E2) pairs for a glass–epoxy ply,
                                  digitised from Tsai (1964), NASA CR-71
```

## Reproducing the manuscript tables

| Table | Script | What it shows |
|---|---|---|
| Table 6 | `run_all` block 1 | RJMCMC vs single-stage MCMC on the experimental data |
| Table S1 | `validation/jcm_revision_runs.m` | effect of the prior on the noise standard deviation |
| Table S2 | `validation/full_validation.m` | 20-replicate study: bias, RMSE, coverage, model identification |
| Table S3 | `validation/jcm_revision_runs.m` | single-stage vs two-stage MCMC and the NLS benchmark |
| Table S4 | `validation/run_model_selection.m` | posterior model probability vs prior on E2f |

Settings used throughout: 10 % burn-in, Gaussian priors (values in the manuscript), jump
probability 0.9 for the model-switching chain, one Metropolis update per parameter per
iteration.

## Data

`data/E2podaci.dat` holds the 15 digitised points (two whitespace-separated columns: fibre
volume fraction `Vf` and transverse modulus `E2` in GPa). The values in the published Table 5
are these numbers rounded to two decimals.

**Provenance.** The points were digitised by the authors from

> Tsai SW. *Structural Behavior of Composite Materials*. NASA CR-71. Washington, DC: National
> Aeronautics and Space Administration; 1964.

which is a US Government work in the public domain. The file is included with this attribution
so that the scripts run out of the box. If you use these numbers, please cite the original
report (Reference 4 in the paper) as well as this repository.

## License

- **Code** (`*.m`): MIT — see `LICENSE`.
- **Digitised data** (`data/E2podaci.dat`): CC0 1.0 Universal (public domain dedication) —
  see `data/NOTICE.txt`.

## Contact

Ivo Marinić-Kragić — imarinic@fesb.hr
Stipe Perišić — sperisic@fesb.hr
University of Split, FESB
