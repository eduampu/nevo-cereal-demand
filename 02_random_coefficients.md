Nevo (2000) replication, part 2: random-coefficients logit
================
Eduardo Ampudia
2026-11-19

## 1. Goal

Part 1 showed that the logit forces substitution to follow market
shares. This notebook estimates the random-coefficients logit of Nevo
(2000), where households differ in how much they value price, sugar and
texture, and compares its substitution patterns with the logit.

Unlike part 1, the estimation here is not coded by hand. It uses the
`BLPestimatoR` package, which ships with the same cereal data. My
contribution in this part is setting up the model, checking the output
and interpreting it.

## 2. What the estimator does

In part 1 a formula turned shares into mean utilities $\delta$, and the
rest was a regression. With household-specific tastes that formula no
longer exists, so the package searches:

1.  **Guess** how much tastes vary across households. These are 13
    numbers: for example, how price sensitivity changes with income, or
    how the taste for sugar changes with age.
2.  **Find $\delta$ by trial and error.** If the model predicts too
    small a share for a brand, raise its $\delta$, and repeat until
    predicted shares match the observed ones.
3.  **Run the regression of part 1** on that $\delta$ (price and brand
    dummies, with the same instruments) and keep the residual $\xi$.
4.  **Score the guess.** If $\xi$ is correlated with the instruments,
    the guess is bad. This score is the GMM objective.
5.  **Improve the guess** and repeat until the score cannot go lower.

The model behind step 1 is

$$u_{ijt} = \delta_{jt} + \mu_{ijt} + \varepsilon_{ijt},$$

where $\mu_{ijt}$ is how household $i$’s valuation of product $j$
differs from the average, given its income, age, children, and an
unobserved taste draw. A household that likes sugar has a high $\mu$ for
every sugary cereal at once, which is what lets similar products be
closer substitutes.

## 3. Data

The package includes Nevo’s cereal data under different column names. As
a check that it is the same dataset as in part 1, I re-run the IV logit
with the function from part 1. The price coefficient should again be
about -30.1.

``` r
dim(productData_cereal)
```

    ## [1] 2256   28

``` r
head(productData_cereal[, c("cdid", "product_id", "share", "price", "sugar", "mushy")])
```

    ##       cdid product_id       share      price sugar mushy
    ## 1 market_1   cereal_1 0.012417212 0.07208794     2     1
    ## 2 market_1   cereal_2 0.007809387 0.11417849    18     1
    ## 3 market_1   cereal_3 0.012994511 0.13239066     4     1
    ## 4 market_1   cereal_4 0.005769961 0.13034408     3     0
    ## 5 market_1   cereal_5 0.017934141 0.15482331    12     0
    ## 6 market_1   cereal_6 0.026601892 0.13704921    14     0

``` r
fit_linear <- function(y, X, Z = X) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  X_hat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(X_hat, X))
  beta  <- bread %*% crossprod(X_hat, y)
  setNames(as.vector(beta), colnames(X))
}

cereal <- productData_cereal
cereal$outside_share <- 1 - ave(cereal$share, cereal$cdid, FUN = sum)
cereal$delta_logit   <- log(cereal$share) - log(cereal$outside_share)

iv_cols       <- paste0("IV", 1:20)
brand_dummies <- model.matrix(~ 0 + factor(productdummy), data = cereal)
logit_iv <- fit_linear(cereal$delta_logit,
                       X = cbind(price = cereal$price, brand_dummies),
                       Z = cbind(as.matrix(cereal[iv_cols]), brand_dummies))
alpha_logit <- -unname(logit_iv["price"])
alpha_logit
```

    ## [1] 30.09775

## 4. Setting up the model

The formula has four parts separated by `|`: variables with a mean
coefficient, exogenous variables, variables with a random coefficient,
and instruments.

``` r
nevos_model <- as.formula(paste(
  "share ~ price + productdummy |",
  "0 + productdummy |",
  "price + sugar + mushy |",
  "0 +", paste(iv_cols, collapse = " + ")
))

# Starting values for delta and the household draws shipped with the package
productData_cereal$startingGuessesDelta <- c(log(w_guesses_cereal))
names(originalDraws_cereal)[1] <- "(Intercept)"

cereal_data <- BLP_data(
  model = nevos_model,
  market_identifier = "cdid",
  par_delta = "startingGuessesDelta",
  product_identifier = "product_id",
  productData = productData_cereal,
  demographic_draws = demographicData_cereal,
  blp_inner_tol = 1e-6,
  blp_inner_maxit = 5000,
  integration_draws = originalDraws_cereal,
  integration_weights = rep(1 / 20, 20)
)
```

The 13 parameters are arranged in a table. Rows are the characteristics
with a random coefficient. The first column is the spread of unobserved
tastes and the others are interactions with demographics. `NA` means the
interaction is left out, as in the paper. These are the starting values
of the search.

``` r
theta_start <- rbind(
  c(0.3302,  5.4819,      NA,  0.2037,     NA),
  c(2.4526, 15.8935, -1.2000,      NA, 2.6342),
  c(0.0163, -0.2506,      NA,  0.0511,     NA),
  c(0.2441,  1.2650,      NA, -0.8091,     NA)
)
rownames(theta_start) <- c("(Intercept)", "price", "sugar", "mushy")
colnames(theta_start) <- c("unobs_sd", "income", "incomesq", "age", "child")
theta_start
```

    ##             unobs_sd  income incomesq     age  child
    ## (Intercept)   0.3302  5.4819       NA  0.2037     NA
    ## price         2.4526 15.8935     -1.2      NA 2.6342
    ## sugar         0.0163 -0.2506       NA  0.0511     NA
    ## mushy         0.2441  1.2650       NA -0.8091     NA

## 5. Estimation

``` r
cereal_est <- estimateBLP(
  blp_data = cereal_data,
  par_theta2 = theta_start,
  solver_method = "BFGS",
  solver_maxit = 1000,
  solver_reltol = 1e-6,
  standardError = "heteroskedastic",
  extremumCheck = FALSE,
  printLevel = 0
)
```

    ## blp_data were prepared with the following arguments:
    ## BLP_data(model = nevos_model, market_identifier = "cdid", product_identifier = "product_id", 
    ##     par_delta = "startingGuessesDelta", productData = productData_cereal, 
    ##     demographic_draws = demographicData_cereal, integration_draws = originalDraws_cereal, 
    ##     integration_weights = rep(1/20, 20), blp_inner_tol = 1e-06, 
    ##     blp_inner_maxit = 5000)
    ## ------------------------------------------ 
    ## Solver message: Successful convergence 
    ## ------------------------------------------ 
    ## Final GMM evaluation at optimal parameters: 
    ## gmm objective: 4.5744
    ##   theta (RC): 0.55 3.29 -0.01 0.09 
    ##   theta (demogr.): 2.3 577.44 -0.38 0.83 0 -29.63 0 0 1.27 0 0.05 -1.35 0 11.03 0 0 
    ##   inner iterations: 70 
    ##   gradient: -0.0605 0.0017 -0.162 -0.0354 -0.0926 -0.0083 -2.2061 0.2945 -0.1468 0.0797 -1.2189 0.2315 2e-04 
    ## Using the heteroskedastic asymptotic variance-covariance matrix...

``` r
summary(cereal_est)
```

    ## 
    ## Data information:
    ## 
    ##   94 market(s) with 2256 products 
    ##   25 linear coefficient(s) (24 exogenous coefficients) 
    ##   13 non-linear parameters related to random coefficients 
    ##   4 demographic variable(s) 
    ## 
    ## Estimation results:
    ## 
    ##  Linear Coefficients
    ##                          Estimate Std. Error   t value     Pr(>|t|)
    ## (Intercept)            -2.5769264  0.8531569 -3.020460 2.523909e-03
    ## price                 -62.1403349 14.5135858 -4.281529 1.856137e-05
    ## productdummyproduct10   1.9000897  0.6653478  2.855784 4.293071e-03
    ## productdummyproduct11   3.6767952  0.6664428  5.517045 3.447468e-08
    ## productdummyproduct12   0.2615993  0.2287542  1.143583 2.527967e-01
    ## productdummyproduct13   2.3591835  0.2184213 10.801071 3.402135e-27
    ## productdummyproduct14   3.8446786  0.6855662  5.608034 2.046376e-08
    ## productdummyproduct15   4.2873240  0.6779154  6.324276 2.544228e-10
    ## productdummyproduct16   6.0916947  0.6723137  9.060792 1.295068e-19
    ## productdummyproduct17   2.6036321  0.6510521  3.999115 6.357976e-05
    ## productdummyproduct18   1.4324858  0.6079004  2.356448 1.845065e-02
    ## productdummyproduct19   4.3135147  0.6191372  6.966977 3.238237e-12
    ## productdummyproduct2    4.8750516  0.6334282  7.696297 1.400664e-14
    ## productdummyproduct20   4.1706988  0.6075700  6.864557 6.669789e-12
    ## productdummyproduct21   4.9106550  0.6827987  7.191951 6.387201e-13
    ## productdummyproduct22   1.7215961  0.5917418  2.909370 3.621574e-03
    ## productdummyproduct23   3.2960677  0.6839398  4.819237 1.441085e-06
    ## productdummyproduct24   1.8879979  0.7307294  2.583717 9.774202e-03
    ## productdummyproduct3    2.3028803  0.2314068  9.951652 2.480314e-23
    ## productdummyproduct4    1.3668373  0.5940217  2.300989 2.139225e-02
    ## 
    ## ...
    ## 
    ##  5 estimates are omitted. They are available in the LinCoefficients generated by summary.
    ## 
    ##  Random Coefficients
    ##                           Estimate   Std. Error    t value     Pr(>|t|)
    ## unobs_sd*(Intercept)   0.551339802   0.16015070  3.4426313 0.0005760842
    ## unobs_sd*price         3.285559241   1.30636137  2.5150462 0.0119016778
    ## unobs_sd*sugar        -0.005237547   0.01336128 -0.3919943 0.6950624576
    ## unobs_sd*mushy         0.091406886   0.18480221  0.4946201 0.6208683096
    ## income*(Intercept)     2.303677738   1.19131427  1.9337280 0.0531465816
    ## income*price         577.439894533 264.92313023  2.1796507 0.0292833612
    ## income*sugar          -0.384035468   0.11965304 -3.2095755 0.0013293115
    ## income*mushy           0.826450256   0.77955170  1.0601609 0.2890713885
    ## incomesq*price       -29.627526224  13.80993187 -2.1453782 0.0319226239
    ## age*(Intercept)        1.268858975   0.62894244  2.0174485 0.0436487313
    ## age*sugar              0.051714014   0.02546546  2.0307513 0.0422802276
    ## age*mushy             -1.350839888   0.66343728 -2.0361230 0.0417380080
    ## child*price           11.025866515   4.11565257  2.6790081 0.0073840611
    ## 
    ##  Wald Test
    ## 126.7789 on  13 DF, p-value: 9.13936382570512e-21 
    ## 
    ## Computational Details: 
    ##   Solver converged with 172 iterations to a minimum at 4.5744 .
    ##   Local minima check: NA
    ##   stopping criterion outer loop: 1e-06
    ##   stopping criterion inner loop: 1e-06
    ##   Market shares are integrated with provided_by_user and 20 draws. 
    ##   Method for standard errors: heteroskedastic

The output has two blocks.

- **Linear coefficients:** the mean price coefficient and the brand
  dummies, the same objects as in part 1. Each brand dummy is that
  brand’s average appeal relative to the brand left out. Five of them
  are cut from the printout.
- **Random coefficients:** the 13 parameters that let households differ.

The last lines say the search converged and at what value of the
objective. Lower is better; on its own the number has no meaning.

### The estimates as a table

``` r
theta_hat <- theta_start
theta_hat[] <- NA
for (i in seq_along(cereal_est$theta_rc)) {
  theta_hat[cereal_est$indices[i, 1], cereal_est$indices[i, 2]] <- cereal_est$theta_rc[i]
}
knitr::kable(round(theta_hat, 3))
```

|             | unobs_sd |  income | incomesq |    age |  child |
|:------------|---------:|--------:|---------:|-------:|-------:|
| (Intercept) |    0.551 |   2.304 |       NA |  1.269 |     NA |
| price       |    3.286 | 577.440 |  -29.628 |     NA | 11.026 |
| sugar       |   -0.005 |  -0.384 |       NA |  0.052 |     NA |
| mushy       |    0.091 |   0.826 |       NA | -1.351 |     NA |

Read each cell as: how the taste for the **row** changes with the
**column**. The demographics are measured relative to the market
average, so a household with average income, age and children has a
taste equal to the mean.

### What each parameter means

``` r
est <- function(name) round(cereal_est$theta_rc[[name]], 3)
mean_price <- round(cereal_est$theta_lin["price", ], 1)

meanings <- data.frame(
  Parameter = c("price (mean)",
                "unobs_sd*price", "income*price", "incomesq*price", "child*price",
                "unobs_sd*(Intercept)", "income*(Intercept)", "age*(Intercept)",
                "unobs_sd*sugar", "income*sugar", "age*sugar",
                "unobs_sd*mushy", "income*mushy", "age*mushy"),
  Estimate = c(mean_price,
               est("unobs_sd*price"), est("income*price"), est("incomesq*price"), est("child*price"),
               est("unobs_sd*(Intercept)"), est("income*(Intercept)"), est("age*(Intercept)"),
               est("unobs_sd*sugar"), est("income*sugar"), est("age*sugar"),
               est("unobs_sd*mushy"), est("income*mushy"), est("age*mushy")),
  Meaning = c(
    "Price sensitivity of a household with average demographics. About twice the logit value of 30.",
    "Spread of price sensitivity that demographics do not explain. Small next to the demographic effects.",
    "With income squared below: how price sensitivity changes with income. Read the two together (see below).",
    "Bends the income effect, so it can be strong at low incomes and fade at high ones.",
    "Positive: holding income fixed, households with children are a little less price sensitive.",
    "Households differ in how much they like cereal in general compared with not buying it.",
    "Positive: richer households are more likely to buy cereal at all.",
    "Positive: older households are more likely to buy cereal at all.",
    "About zero: no taste differences for sugar beyond those explained by income and age.",
    "Negative: richer households like sugary cereals less.",
    "Positive: older households like sugar more, as measured here.",
    "Small and imprecise: little unexplained variation in taste for soggy cereals.",
    "Positive but imprecise: weak evidence that richer households like soggy cereals more.",
    "Negative: older households dislike cereals that go soggy in milk."
  )
)
knitr::kable(meanings, row.names = FALSE)
```

| Parameter | Estimate | Meaning |
|:---|---:|:---|
| price (mean) | -62.100 | Price sensitivity of a household with average demographics. About twice the logit value of 30. |
| unobs_sd\*price | 3.286 | Spread of price sensitivity that demographics do not explain. Small next to the demographic effects. |
| income\*price | 577.440 | With income squared below: how price sensitivity changes with income. Read the two together (see below). |
| incomesq\*price | -29.628 | Bends the income effect, so it can be strong at low incomes and fade at high ones. |
| child\*price | 11.026 | Positive: holding income fixed, households with children are a little less price sensitive. |
| unobs_sd\*(Intercept) | 0.551 | Households differ in how much they like cereal in general compared with not buying it. |
| income\*(Intercept) | 2.304 | Positive: richer households are more likely to buy cereal at all. |
| age\*(Intercept) | 1.269 | Positive: older households are more likely to buy cereal at all. |
| unobs_sd\*sugar | -0.005 | About zero: no taste differences for sugar beyond those explained by income and age. |
| income\*sugar | -0.384 | Negative: richer households like sugary cereals less. |
| age\*sugar | 0.052 | Positive: older households like sugar more, as measured here. |
| unobs_sd\*mushy | 0.091 | Small and imprecise: little unexplained variation in taste for soggy cereals. |
| income\*mushy | 0.826 | Positive but imprecise: weak evidence that richer households like soggy cereals more. |
| age\*mushy | -1.351 | Negative: older households dislike cereals that go soggy in milk. |

A parameter is estimated precisely when its estimate is large relative
to its standard error (a t value above about 2 in the summary).
`unobs_sd*sugar`, `unobs_sd*mushy` and `income*mushy` are not, so their
signs should not be read much into.

### Price sensitivity across households

The income parameters are hard to read on their own: in this data income
and income squared are almost perfectly correlated, so a large positive
coefficient on one and a large negative one on the other partly cancel.
The clearer way to read them is to compute each simulated household’s
own price coefficient and look at the result.

``` r
draws <- function(x) as.matrix(x[, grep("^draw_", names(x))])

household_price_coef <- cereal_est$theta_lin["price", ] +
  theta_hat["price", "unobs_sd"] * draws(originalDraws_cereal$price) +
  theta_hat["price", "income"]   * draws(demographicData_cereal$income) +
  theta_hat["price", "incomesq"] * draws(demographicData_cereal$incomesq) +
  theta_hat["price", "child"]    * draws(demographicData_cereal$child)

income <- draws(demographicData_cereal$income)
income_group <- cut(income, quantile(income, 0:4 / 4), include.lowest = TRUE,
                    labels = c("Lowest 25%", "Second 25%", "Third 25%", "Highest 25%"))

knitr::kable(data.frame(
  Statistic = c("Median household", "10th percentile (most price sensitive)",
                "90th percentile (least price sensitive)", "Mean",
                "Share with a positive coefficient"),
  Value = c(median(household_price_coef), quantile(household_price_coef, 0.1),
            quantile(household_price_coef, 0.9), mean(household_price_coef),
            mean(household_price_coef > 0))
), digits = 2)
```

| Statistic                               |   Value |
|:----------------------------------------|--------:|
| Median household                        |  -34.48 |
| 10th percentile (most price sensitive)  | -126.35 |
| 90th percentile (least price sensitive) |  -21.86 |
| Mean                                    |  -62.20 |
| Share with a positive coefficient       |    0.00 |

``` r
knitr::kable(data.frame(
  `Income group` = levels(income_group),
  `Median price coefficient` = as.vector(tapply(household_price_coef, income_group, median)),
  check.names = FALSE
), digits = 1)
```

| Income group | Median price coefficient |
|:-------------|-------------------------:|
| Lowest 25%   |                   -104.9 |
| Second 25%   |                    -39.9 |
| Third 25%    |                    -26.9 |
| Highest 25%  |                    -27.4 |

**What these numbers say.**

- The median household has a price coefficient of -34.5, close to the
  logit’s 30. The mean of -62.2 is pulled down by a minority of very
  price-sensitive households.
- Those households are the poorest: the lowest income quarter has a
  median of -104.9, while the top half is around -27.4. Price
  sensitivity falls with income and then flattens, which is the curved
  shape that income squared allows.
- No household has a positive coefficient, so demand slopes down for
  everyone.

### Comparison with the paper

These results match the PyBLP replication of the same model on the same
data, which reports an objective of about 4.6, a mean price coefficient
of -62.7, and interactions of about 588 with income and -30.2 with
income squared. PyBLP’s authors describe these results as similar to
those in the original paper. The small gaps with the numbers here (for
example 577 against 588) come from the optimiser stopping slightly
before the exact minimum.

## 6. Elasticities

The package computes elasticities from the estimates. Rows are the
product whose share responds and columns the product whose price
changes, as in part 1.

``` r
price_coef <- cereal_est$theta_lin["price", ]

# Store the estimated mean utilities in the data object
delta_data <- data.frame(
  product_id = cereal_data$parameters$product_id,
  cdid = cereal_data$parameters$market_id_char_in,
  startingGuessesDelta = cereal_est$delta
)
cereal_data_updated <- update_BLP_data(data_update = delta_data, blp_data = cereal_data)
```

    ## Mean utility variable startingGuessesDelta has been updated.

``` r
share_info <- getShareInfo(blp_data = cereal_data_updated, par_theta2 = theta_hat,
                           printLevel = 0)

market_elasticities <- function(market) {
  get_elasticities(blp_data = cereal_data_updated, share_info = share_info,
                   theta_lin = price_coef, variable = "price",
                   market = market, printLevel = 0)
}

# Prices and shares of a market, in the same product order as the elasticity matrix
market_products <- function(market, E) {
  m <- cereal[cereal$cdid == market, ]
  if (!is.null(rownames(E)) && all(rownames(E) %in% m$product_id)) {
    m <- m[match(rownames(E), m$product_id), ]
  }
  m
}
```

### Same table as part 1

``` r
E_rc <- market_elasticities("market_1")
m1   <- market_products("market_1", E_rc)

E_logit <- matrix(alpha_logit * m1$price * m1$share,
                  nrow = nrow(m1), ncol = nrow(m1), byrow = TRUE)
diag(E_logit) <- -alpha_logit * m1$price * (1 - m1$share)
dimnames(E_logit) <- list(m1$product_id, m1$product_id)

knitr::kable(E_rc[1:5, 1:5], digits = 3, caption = "Random coefficients")
```

|          | cereal_1 | cereal_2 | cereal_3 | cereal_4 | cereal_5 |
|:---------|---------:|---------:|---------:|---------:|---------:|
| cereal_1 |   -2.337 |    0.008 |    0.130 |    0.055 |    0.041 |
| cereal_2 |    0.008 |   -4.615 |    0.028 |    0.013 |    0.136 |
| cereal_3 |    0.067 |    0.015 |   -3.596 |    0.042 |    0.068 |
| cereal_4 |    0.066 |    0.016 |    0.095 |   -3.999 |    0.094 |
| cereal_5 |    0.013 |    0.044 |    0.042 |    0.025 |   -4.974 |

Random coefficients

``` r
knitr::kable(E_logit[1:5, 1:5], digits = 3, caption = "Logit")
```

|          | cereal_1 | cereal_2 | cereal_3 | cereal_4 | cereal_5 |
|:---------|---------:|---------:|---------:|---------:|---------:|
| cereal_1 |   -2.143 |    0.027 |    0.052 |    0.023 |    0.084 |
| cereal_2 |    0.027 |   -3.410 |    0.052 |    0.023 |    0.084 |
| cereal_3 |    0.027 |    0.027 |   -3.933 |    0.023 |    0.084 |
| cereal_4 |    0.027 |    0.027 |    0.052 |   -3.900 |    0.084 |
| cereal_5 |    0.027 |    0.027 |    0.052 |    0.023 |   -4.576 |

Logit

``` r
knitr::kable(m1[1:5, c("product_id", "sugar", "mushy", "price", "share")],
             digits = 3, row.names = FALSE, caption = "Characteristics of these products")
```

| product_id | sugar | mushy | price | share |
|:-----------|------:|------:|------:|------:|
| cereal_1   |     2 |     1 | 0.072 | 0.012 |
| cereal_2   |    18 |     1 | 0.114 | 0.008 |
| cereal_3   |     4 |     1 | 0.132 | 0.013 |
| cereal_4   |     3 |     0 | 0.130 | 0.006 |
| cereal_5   |    12 |     0 | 0.155 | 0.018 |

Characteristics of these products

In the logit table every column repeats one number. With random
coefficients the numbers differ within a column: how much a rival gains
depends on which rival it is, not only on the product raising its price.

**An example from the table.** `cereal_3` has little sugar and goes
soggy in milk. When its price rises 1%, `cereal_1`, which is similar
(low sugar, soggy), gains 0.13% of sales, while `cereal_2`, which is
very sugary, gains only 0.028%. In the logit both would gain 0.052%. The
model has learned that similar cereals compete more closely.

The own-price elasticities on the diagonal are similar in the two
models. The average sensitivity of each brand to its own price barely
changes; what changes is where the lost customers go.

### Diversion ratios

The diversion ratio from $k$ to $j$ is the fraction of the customers $k$
loses after a price increase that go to $j$. It is the first statistic
looked at in a merger review: if two merging products divert a lot to
each other, the merged firm recovers sales it would otherwise lose, so
it can raise prices.

An elasticity is a percentage change. To count customers, I convert it
back to a change in share:
$\partial s_j / \partial p_k = \eta_{jk} \, s_j / p_k$.

``` r
diversion <- function(E, s, p) {
  derivs <- E * outer(s, 1 / p)            # change in share j per unit of price k
  out <- -sweep(derivs, 2, diag(derivs), "/")
  diag(out) <- 0
  out                                      # [j, k]: share of k's lost customers that go to j
}

Div_rc    <- diversion(E_rc, m1$share, m1$price)
Div_logit <- diversion(E_logit, m1$share, m1$price)

knitr::kable(Div_rc[1:5, 1:5], digits = 3,
             caption = "Random coefficients: diversion from column to row")
```

|          | cereal_1 | cereal_2 | cereal_3 | cereal_4 | cereal_5 |
|:---------|---------:|---------:|---------:|---------:|---------:|
| cereal_1 |    0.000 |    0.003 |    0.034 |    0.030 |    0.006 |
| cereal_2 |    0.002 |    0.000 |    0.005 |    0.004 |    0.012 |
| cereal_3 |    0.030 |    0.005 |    0.000 |    0.024 |    0.010 |
| cereal_4 |    0.013 |    0.003 |    0.012 |    0.000 |    0.006 |
| cereal_5 |    0.008 |    0.022 |    0.016 |    0.020 |    0.000 |

Random coefficients: diversion from column to row

``` r
knitr::kable(Div_logit[1:5, 1:5], digits = 3,
             caption = "Logit: diversion from column to row")
```

|          | cereal_1 | cereal_2 | cereal_3 | cereal_4 | cereal_5 |
|:---------|---------:|---------:|---------:|---------:|---------:|
| cereal_1 |    0.000 |    0.013 |    0.013 |    0.012 |    0.013 |
| cereal_2 |    0.008 |    0.000 |    0.008 |    0.008 |    0.008 |
| cereal_3 |    0.013 |    0.013 |    0.000 |    0.013 |    0.013 |
| cereal_4 |    0.006 |    0.006 |    0.006 |    0.000 |    0.006 |
| cereal_5 |    0.018 |    0.018 |    0.018 |    0.018 |    0.000 |

Logit: diversion from column to row

**How to read one number.** Each column is a brand that raises its
price, and the rows are where its lost customers go. Under random
coefficients, 3% of the customers `cereal_1` loses go to `cereal_3`,
against 1.3% in the logit. In the logit every number in a row is the
same, because diversion depends only on the size of the brand receiving
customers.

Whatever does not go to another brand goes to the outside good, that is,
the customer stops buying cereal. Averaging over all products and
markets:

``` r
markets <- unique(cereal$cdid)
by_market <- lapply(markets, function(market) {
  E <- market_elasticities(market)
  m <- market_products(market, E)
  data.frame(own_rc     = diag(E),
             own_logit  = -alpha_logit * m$price * (1 - m$share),
             outside_rc = 1 - colSums(diversion(E, m$share, m$price)),
             outside_logit = m$outside_share / (1 - m$share))
})
by_market <- do.call(rbind, by_market)

comparison <- data.frame(
  Model = c("Logit", "Random coefficients"),
  `Mean own-price elasticity` = c(mean(by_market$own_logit), mean(by_market$own_rc)),
  `Share of products with inelastic demand` =
    c(mean(by_market$own_logit > -1), mean(by_market$own_rc > -1)),
  `Mean diversion to the outside good` =
    c(mean(by_market$outside_logit), mean(by_market$outside_rc)),
  check.names = FALSE
)
knitr::kable(comparison, digits = 3, row.names = FALSE)
```

| Model | Mean own-price elasticity | Share of products with inelastic demand | Mean diversion to the outside good |
|:---|---:|---:|---:|
| Logit | -3.713 | 0 | 0.535 |
| Random coefficients | -3.621 | 0 | 0.367 |

**What these numbers say.**

- **Own-price elasticity:** about -3.6 in both models. On average, a 1%
  price rise loses a brand a bit under 4% of its sales either way.
- **Inelastic products:** none in either model, as profit-maximising
  pricing requires.
- **Diversion to the outside good:** 53% of lost customers stop buying
  cereal in the logit, against 37% with random coefficients. The logit
  sends too many customers out of the market, because the outside good
  is the biggest “product” and the logit spreads customers by size.

**Why does this matter?** In the logit, a large share of the customers
who leave a brand stop buying cereal. With random coefficients more of
them stay in the market and move to another brand. Brands compete with
each other more closely than the logit implies, so a merger between
close substitutes has a larger effect on prices.

## 7. Saved for part 3

Part 3 needs the derivatives of every share with respect to every price
under this model. They are saved here so part 3 does not have to repeat
the estimation.

``` r
dir.create("output", showWarnings = FALSE)
rc_derivatives <- lapply(markets, function(market) {
  E <- market_elasticities(market)
  m <- market_products(market, E)
  list(market = market, product_id = m$product_id, price = m$price, share = m$share,
       derivs = E * outer(m$share, 1 / m$price))   # [j, k] = d s_j / d p_k
})
saveRDS(rc_derivatives, file.path("output", "rc_derivatives.rds"))
```

## References

- Berry, S., Levinsohn, J., and Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica* 63(4).
- Brunner, D. (maintainer). BLPestimatoR: Performs a BLP Demand
  Estimation. R package, CRAN.
- Conlon, C., and Gortmaker, J. (2020). Best Practices for
  Differentiated Products Demand Estimation with PyBLP. *RAND Journal of
  Economics* 51(4).
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy* 9(4).
