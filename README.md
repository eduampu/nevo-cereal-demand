# Demand estimation for differentiated products: replicating Nevo (2000)

Replication in R of the numerical example in Nevo, A. (2000), "A Practitioner's
Guide to Estimation of Random-Coefficients Logit Models of Demand", *Journal of
Economics & Management Strategy* 9(4), 513-548, applied to the ready-to-eat
cereal market. The question behind it is the one a merger review asks: when a
brand raises its price, where do its customers go?

## Contents

| File | What it does | Status |
|---|---|---|
| `01_logit.Rmd` | Logit demand estimated by hand with matrix algebra (OLS and IV with brand fixed effects), checked against `fixest`. Elasticities and why the logit is not enough for merger analysis. | Done |
| `02_random_coefficients.Rmd` | Random-coefficients logit estimated with the `BLPestimatoR` package. Comparison with the PyBLP replication, elasticities and diversion ratios against the logit. | Draft |
| `03_merger_simulation.Rmd` | Simulated merger of the two largest firms: marginal costs from profit maximisation, post-merger prices and consumer losses under the logit, and pricing pressure (GUPPI) under both demand models. | Draft |

The `.md` files are the knitted notebooks and can be read directly on GitHub.

## Main results

**Demand (parts 1 and 2)**

- Controlling for brand quality and instrumenting price changes the logit price
  coefficient from -10.1 (OLS) to -30.1 (IV with brand fixed effects). A 1 cent
  increase in the price of a serving cuts a small brand's sales by about a
  quarter.
- The average own-price elasticity is about -3.7 in the logit and -3.6 with
  random coefficients, and no product has inelastic demand in either model.
- The random-coefficients estimates match the PyBLP replication of Nevo's model
  (mean price coefficient about -62). The median household is about as price
  sensitive as in the logit; the lower mean comes from the poorest quarter of
  households, which is far more price sensitive than the rest.
- The models differ in where customers go. When a brand raises its price, the
  logit says 54% of the customers it loses stop buying cereal; random
  coefficients say 37%, with more of them moving to similar brands.

**Merger of the two largest firms (part 3)**

- Before the merger, the median price is about 31% above marginal cost.
- With logit demand, the merging firms raise prices by a median of 5.4% and
  6.3%, rivals by under 0.1%, and consumers lose the equivalent of about 4.2%
  of their cereal spending.
- The pricing-pressure index (GUPPI) is 5.6% under the logit and 8.7% under
  random coefficients, because the richer model finds that more of each merging
  brand's lost customers go to its partner (24% against 17%). The logit
  simulation is therefore likely to understate the price effect.

## Data

The data are the simulated ("fake") cereal data that Nevo distributed with the
paper; the original scanner data are proprietary.

- 94 markets (47 cities x 2 quarters), 24 products per market, 2,256 observations.
- 20 simulated households per market with demographics (income, age, child).
- Part 1 downloads the files from the [PyBLP](https://github.com/jeffgortmaker/pyblp)
  repository (Conlon and Gortmaker, 2020). Part 2 uses the copy included in
  `BLPestimatoR`, and checks that both give the same logit estimate.

## How to run

```r
install.packages(c("rmarkdown", "dplyr", "knitr", "fixest", "BLPestimatoR"))
rmarkdown::render("01_logit.Rmd")
rmarkdown::render("02_random_coefficients.Rmd")   # also saves output/ used by part 3
rmarkdown::render("03_merger_simulation.Rmd")
```

## How this was made

I built this project with the help of an AI assistant (Claude), which wrote
most of the code under my direction. I ran everything, checked the results
(part 1 against `fixest`, part 2 against the logit and the PyBLP
replication), and worked through the economics until I could explain each step.
The notes in the notebooks are that explanation.

## References

- Berry, S. (1994). Estimating Discrete-Choice Models of Product Differentiation. *RAND Journal of Economics* 25(2).
- Berry, S., Levinsohn, J., and Pakes, A. (1995). Automobile Prices in Market Equilibrium. *Econometrica* 63(4).
- Conlon, C., and Gortmaker, J. (2020). Best Practices for Differentiated Products Demand Estimation with PyBLP. *RAND Journal of Economics* 51(4).
- Nevo, A. (2000). A Practitioner's Guide to Estimation of Random-Coefficients Logit Models of Demand. *Journal of Economics & Management Strategy* 9(4).
