Nevo (2000) replication, part 1: logit demand
================
Eduardo Ampudia
2026-11-19

## 1. Goal

This notebook estimates the plain logit demand model on Nevo’s cereal
data. It is the first step of the replication and the benchmark that the
random-coefficients model in part 2 improves on. The estimator is
computed with matrix algebra and then checked against `fixest`.

## 2. Data

The files come from the PyBLP repository. They are downloaded once and
kept in `data/`.

``` r
data_dir <- "data"
base_url <- "https://raw.githubusercontent.com/jeffgortmaker/pyblp/master/pyblp/data/"
dir.create(data_dir, showWarnings = FALSE)
for (f in c("nevo_products.csv", "nevo_agents.csv")) {
  dest <- file.path(data_dir, f)
  if (!file.exists(dest)) download.file(paste0(base_url, f), dest, mode = "wb")
}
products <- read.csv(file.path(data_dir, "nevo_products.csv"), stringsAsFactors = FALSE)
```

``` r
required <- c("market_ids", "product_ids", "firm_ids", "shares", "prices", "sugar", "mushy")
not_found <- setdiff(required, names(products))
if (length(not_found) > 0) {
  stop("Columns not found: ", paste(not_found, collapse = ", "),
       ". Columns in the file: ", paste(names(products), collapse = ", "))
}
iv_cols <- grep("^demand_instruments", names(products), value = TRUE)

data_summary <- data.frame(
  Observations = nrow(products),
  Markets      = n_distinct(products$market_ids),
  Products     = n_distinct(products$product_ids),
  Firms        = n_distinct(products$firm_ids),
  Instruments  = length(iv_cols)
)
knitr::kable(data_summary)
```

| Observations | Markets | Products | Firms | Instruments |
|-------------:|--------:|---------:|------:|------------:|
|         2256 |      94 |       24 |     5 |          20 |

**What these numbers say.** 2,256 rows = 24 brands × 94 markets, so
every brand is observed in every market. A market is one city in one
quarter. The 24 brands belong to 5 firms, which matters in part 3. The
20 instruments are used in section 4.

Each row is a cereal brand in a city-quarter. `shares` is the brand’s
share of the potential market (one serving per person per day), so
shares in a market add up to less than one. The remainder is the outside
good: not buying any of the 24 brands.

``` r
products <- products %>%
  group_by(market_ids) %>%
  mutate(outside_share = 1 - sum(shares)) %>%
  ungroup() %>%
  mutate(delta = log(shares) - log(outside_share))

summary(products[c("shares", "outside_share", "prices", "sugar", "mushy")])
```

    ##      shares          outside_share        prices            sugar       
    ##  Min.   :0.0001818   Min.   :0.3046   Min.   :0.04549   Min.   : 0.000  
    ##  1st Qu.:0.0051832   1st Qu.:0.4395   1st Qu.:0.10549   1st Qu.: 3.000  
    ##  Median :0.0111408   Median :0.5358   Median :0.12383   Median : 8.500  
    ##  Mean   :0.0198251   Mean   :0.5242   Mean   :0.12574   Mean   : 8.625  
    ##  3rd Qu.:0.0246458   3rd Qu.:0.6043   3rd Qu.:0.14329   3rd Qu.:13.250  
    ##  Max.   :0.4468826   Max.   :0.8152   Max.   :0.22573   Max.   :20.000  
    ##      mushy       
    ##  Min.   :0.0000  
    ##  1st Qu.:0.0000  
    ##  Median :0.0000  
    ##  Mean   :0.3333  
    ##  3rd Qu.:1.0000  
    ##  Max.   :1.0000

**What these numbers say.**

- **shares:** the median brand sells 1.1% of all the servings people in
  its city could eat. Each brand is small on its own.
- **outside_share:** in the median market, 54% of potential servings are
  not one of these 24 brands. The outside good is by far the largest
  “product”, which matters later.
- **prices:** prices are in dollars per serving. The median is \$0.124,
  about 12 cents.
- **sugar and mushy:** sugar content ranges from 0 to
  20. `mushy` is 1 if the cereal goes soggy in milk; a third of the
      brands do.

## 3. Model

Consumer $i$ in market $t$ gets utility from product $j$

$$u_{ijt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt} + \varepsilon_{ijt},$$

where $\xi_{jt}$ is product quality that consumers and firms observe but
the econometrician does not. With extreme-value $\varepsilon_{ijt}$,
market shares have a closed form that can be inverted (Berry, 1994):

$$\ln s_{jt} - \ln s_{0t} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}.$$

The left-hand side is data (`delta` above), so this is a linear
regression with $\xi_{jt}$ as the error term. Prices are endogenous:
firms charge more for products with high $\xi_{jt}$, which biases the
OLS price coefficient towards zero. Nevo addresses this in two ways:

- **Brand fixed effects** absorb the mean quality of each brand, so the
  error is only the market-specific deviation from it.
- **Instruments**: prices of the same brand in other cities, which share
  cost shocks but are assumed not to share demand shocks. These are the
  20 `demand_instruments` columns.

### Notes on the logic

**Why is a market share a probability?** I cannot see
$\varepsilon_{ijt}$, so I cannot say what a given consumer buys, only
how likely each product is to give them the highest utility. With many
consumers, that probability is the share of people who buy the product.

**What is $\delta_{jt}$?** Write utility as
$\delta_{jt} + \varepsilon_{ijt}$, with
$\delta_{jt} = x_{jt}\beta - \alpha p_{jt} + \xi_{jt}$. It is the
average valuation of product $j$ across consumers. It is the same for
everyone but different for each product, so it does not cancel when
products are compared. A higher $\delta_{jt}$ tilts the odds towards
product $j$, and $\varepsilon_{ijt}$ is the reason the best product on
average does not get the whole market. If every product had the same
$\delta$, all shares would be equal.

**Do consumers really agree on price and quality?** No. The logit
assumes they do: $\alpha$ and $\beta$ carry no $i$ subscript, so every
household is equally price sensitive and values sugar equally. Consumers
differ only through $\varepsilon_{ijt}$, which is drawn independently
for each product. It can capture “this person likes this brand” but not
“this person likes sugary cereals”.

**Why the extreme-value distribution?** For convenience, not because
tastes follow it. It is the distribution for which the probability of
being the best option has a closed form,
$s_{jt} = \exp(\delta_{jt}) / (1 + \sum_k \exp(\delta_{kt}))$. Dividing
by the outside share, $s_{0t} = 1 / (1 + \sum_k \exp(\delta_{kt}))$,
cancels the denominator and gives the equation above.

**What does that assumption cost?** The shape of the distribution
matters little. The cost comes from the draws being independent across
products and identical across consumers:

- Substitution depends on market shares, not on how similar two products
  are.
- Own-price elasticities are close to $-\alpha p_{jt}$, so expensive
  products look more elastic by construction.
- Implied markups follow from those elasticities, so they are imposed
  too.

## 4. Estimation by hand

One function covers OLS and two-stage least squares, since OLS is the
case where the instruments are the regressors themselves. Standard
errors are heteroskedasticity-robust (HC0).

``` r
fit_linear <- function(y, X, Z = X) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  X_hat <- Z %*% solve(crossprod(Z), crossprod(Z, X))  # first stage fitted values
  bread <- solve(crossprod(X_hat, X))
  beta  <- bread %*% crossprod(X_hat, y)
  resid <- as.vector(y - X %*% beta)                   # residuals use actual X
  meat  <- crossprod(X_hat * resid)
  V     <- bread %*% meat %*% t(bread)
  list(coef  = setNames(as.vector(beta), colnames(X)),
       se    = setNames(sqrt(diag(V)), colnames(X)),
       resid = resid)
}
```

Three specifications:

``` r
y <- products$delta
brand_dummies <- model.matrix(~ 0 + factor(product_ids), data = products)
instruments   <- as.matrix(products[iv_cols])

X_chars <- cbind(const  = 1,
                 prices = products$prices,
                 sugar  = products$sugar,
                 mushy  = products$mushy)
X_fe <- cbind(prices = products$prices, brand_dummies)
Z_fe <- cbind(instruments, brand_dummies)

ols_chars <- fit_linear(y, X_chars)
ols_fe    <- fit_linear(y, X_fe)
iv_fe     <- fit_linear(y, X_fe, Z_fe)
```

``` r
results <- data.frame(
  Specification   = c("OLS, characteristics", "OLS, brand FE", "IV, brand FE"),
  `Price coef.`   = c(ols_chars$coef["prices"], ols_fe$coef["prices"], iv_fe$coef["prices"]),
  `Robust SE`     = c(ols_chars$se["prices"], ols_fe$se["prices"], iv_fe$se["prices"]),
  `Brand FE`      = c("No", "Yes", "Yes"),
  Instruments     = c("No", "No", "Yes"),
  check.names = FALSE
)
knitr::kable(results, digits = 3, row.names = FALSE)
```

| Specification        | Price coef. | Robust SE | Brand FE | Instruments |
|:---------------------|------------:|----------:|:---------|:------------|
| OLS, characteristics |     -10.120 |     0.831 | No       | No          |
| OLS, brand FE        |     -28.950 |     0.977 | Yes      | No          |
| IV, brand FE         |     -30.098 |     1.019 | Yes      | Yes         |

**What the price coefficient means.** The coefficient is $-\alpha$, the
utility lost per dollar of price. Prices are per serving, so a 1 cent
increase (\$0.01) lowers a brand’s mean utility by
$\alpha \times 0.01 =$ 0.3. For a small brand, that cuts its sales by
about 26%. A 1 cent rise is about an 8% price increase, so cereal buyers
are quite price sensitive.

**Why the three rows differ.**

- **OLS with characteristics (-10.1):** far too small. Expensive brands
  are often the ones people like most, for reasons the data do not
  record (advertising, brand image). OLS sees expensive brands selling
  well and concludes that price barely matters.
- **OLS with brand fixed effects (-28.9):** the brand dummies absorb
  each brand’s average appeal, which removes most of that bias. The
  coefficient almost triples.
- **IV with brand fixed effects (-30.1):** the instruments remove what
  is left, local demand shocks that also move local prices. The change
  is small, so after the brand dummies little bias remains.

The robust standard errors are around 1, so the price coefficient is
estimated precisely.

### Instrument strength

F statistic for the 20 excluded instruments in the first stage,
conditional on the brand fixed effects.

``` r
rss <- function(y, X) {
  X <- as.matrix(X)
  sum((y - X %*% solve(crossprod(X), crossprod(X, y)))^2)
}
rss_restricted   <- rss(products$prices, brand_dummies)
rss_unrestricted <- rss(products$prices, Z_fe)
q <- length(iv_cols)
first_stage_F <- ((rss_restricted - rss_unrestricted) / q) /
  (rss_unrestricted / (nrow(products) - ncol(Z_fe)))
first_stage_F
```

    ## [1] 3363.84

**What this number means.** The F statistic asks whether the instruments
help predict price once the brand is known. The usual rule of thumb for
weak instruments is F below 10. At 3364, the instruments are very
strong: a brand’s price in other cities predicts its local price well.
Strength is not validity, though. The instruments are valid only if
demand shocks are not shared across cities (national advertising would
break this).

### Check against `fixest`

The point estimate has to match exactly. `fixest` applies a small-sample
correction to robust standard errors, so those differ slightly.

``` r
if (requireNamespace("fixest", quietly = TRUE)) {
  iv_formula <- as.formula(paste(
    "delta ~ 1 | product_ids | prices ~", paste(iv_cols, collapse = " + ")
  ))
  fixest_iv <- fixest::feols(iv_formula, data = products, vcov = "hetero")
  comparison <- data.frame(
    Source        = c("By hand", "fixest"),
    `Price coef.` = c(iv_fe$coef["prices"], coef(fixest_iv)["fit_prices"]),
    `Robust SE`   = c(iv_fe$se["prices"], fixest::se(fixest_iv)["fit_prices"]),
    check.names = FALSE
  )
  knitr::kable(comparison, digits = 4, row.names = FALSE)
} else {
  message("fixest is not installed; skipping the check.")
}
```

| Source  | Price coef. | Robust SE |
|:--------|------------:|----------:|
| By hand |    -30.0978 |    1.0187 |
| fixest  |    -30.0978 |    1.0244 |

**What this check shows.** The two price coefficients are identical, so
the matrix code in section 4 is correct. The standard errors differ only
in the third decimal because `fixest` applies a small-sample correction.

## 5. Elasticities

In the logit model the elasticity of product $j$’s share with respect to
the price of product $k$ is

$$\eta_{jk} = \begin{cases} -\alpha\, p_j (1 - s_j) & j = k \\ \alpha\, p_k s_k & j \neq k. \end{cases}$$

``` r
alpha <- -unname(iv_fe$coef["prices"])

logit_elasticities <- function(p, s, alpha) {
  J <- length(p)
  E <- matrix(alpha * p * s, nrow = J, ncol = J, byrow = TRUE)  # E[j, k] = alpha * p_k * s_k
  diag(E) <- -alpha * p * (1 - s)
  E
}

products <- products %>% mutate(own_elasticity = -alpha * prices * (1 - shares))
own_summary <- data.frame(
  Mean   = mean(products$own_elasticity),
  Median = median(products$own_elasticity),
  Min    = min(products$own_elasticity),
  Max    = max(products$own_elasticity),
  `Share inelastic` = mean(products$own_elasticity > -1),
  check.names = FALSE
)
knitr::kable(own_summary, digits = 3)
```

|   Mean | Median |    Min |    Max | Share inelastic |
|-------:|-------:|-------:|-------:|----------------:|
| -3.713 | -3.655 | -6.634 | -1.334 |               0 |

**What these numbers say.**

- **Mean -3.71:** if a typical brand raises its price by 1%, it loses
  about 3.7% of its sales.
- **Range:** from -1.33 to -6.63. In the logit this is almost entirely
  $-\alpha p$, so expensive brands look more price sensitive by
  construction, not because the data say so.
- **Share inelastic = 0:** no brand has an elasticity between 0 and -1.
  That is a good sign. A firm facing inelastic demand could raise its
  price, lose few customers and earn more, so a profit-maximising firm
  never prices there.
- **Link to margins:** for a firm selling one product, the
  profit-maximising margin is $(p - c)/p = 1/|\eta|$. An elasticity of
  -3.7 implies a margin of about 27%. Part 3 uses this logic to recover
  costs.

### Why the logit is not enough for merger analysis

Elasticities in the first market, first five products. Rows are the
product whose share responds, columns the product whose price changes.

``` r
first_market <- products %>% filter(market_ids == market_ids[1])
E <- logit_elasticities(first_market$prices, first_market$shares, alpha)
dimnames(E) <- list(first_market$product_ids, first_market$product_ids)
knitr::kable(E[1:5, 1:5], digits = 3)
```

|       |  F1B04 |  F1B06 |  F1B07 |  F1B09 |  F1B11 |
|:------|-------:|-------:|-------:|-------:|-------:|
| F1B04 | -2.143 |  0.027 |  0.052 |  0.023 |  0.084 |
| F1B06 |  0.027 | -3.410 |  0.052 |  0.023 |  0.084 |
| F1B07 |  0.027 |  0.027 | -3.933 |  0.023 |  0.084 |
| F1B09 |  0.027 |  0.027 |  0.052 | -3.900 |  0.084 |
| F1B11 |  0.027 |  0.027 |  0.052 |  0.023 | -4.576 |

Off the diagonal, every entry in a column is the same number. When a
product raises its price, all rivals gain share in the same proportion,
whether or not they resemble it. A sugary children’s cereal is treated
as an equally close substitute for another children’s cereal as for a
bran cereal. The price effect of a merger depends on how closely the
merging products compete, so this restriction decides the answer before
looking at the data. Part 2 removes it by letting tastes for price,
sugar, and texture vary across consumers.

Two more things to read from the table:

- The cross elasticities are small compared with the own elasticities.
  In the logit, customers who leave a product spread out by share, and
  the outside good has the largest share, so the model sends most of
  them out of the cereal market instead of to a rival brand. That makes
  any merger look harmless.
- The five products shown belong to the same firm (the product code
  starts with the firm number). The model does not treat them as closer
  substitutes for each other than for any other firm’s products.

### Why estimate the logit at all

- **It is a building block.** The random-coefficients estimator recovers
  $\delta_{jt}$ for each guess of the taste parameters and then runs
  this same IV regression inside its search.
- **It gives starting values** for the numerical search in part 2.
- **It isolates the endogeneity problem.** The OLS and IV comparison
  tests the instruments on their own, before adding taste heterogeneity.
- **It is the benchmark.** Simulating the same merger under both models
  shows whether the richer model changes the answer.

## 6. Next

- `02_random_coefficients.Rmd`: random-coefficients logit estimated with
  the `BLPestimatoR` package, compared with the results in the paper and
  with the logit.
- `03_merger_simulation.Rmd`: marginal costs from the firms’ first-order
  conditions and simulated post-merger prices.
