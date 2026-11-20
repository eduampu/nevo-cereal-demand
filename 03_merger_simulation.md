Nevo (2000) replication, part 3: merger simulation
================
Eduardo Ampudia
2026-11-19

## 1. Goal

Parts 1 and 2 estimated how cereal demand responds to prices. This
notebook uses those estimates for the question a competition authority
would ask: if the two largest cereal makers merged, how much would
prices rise?

The method is the standard merger simulation used in antitrust work,
applied to cereal by Nevo in a companion paper (Nevo, 2000, *RAND
Journal of Economics*). The data are Nevo’s simulated data, where firm
and brand names are anonymised, so the merger is hypothetical.

## 2. The theory

### Firms set prices to maximise profit

Each firm chooses the prices of all its brands at once, taking the
prices of other firms as given (Bertrand competition). A firm that owns
brands A and B earns

$$\pi = (p_A - c_A)\, s_A(p) + (p_B - c_B)\, s_B(p),$$

where $c$ is marginal cost and $s(p)$ comes from the estimated demand.
Raising the price of A has three effects:

1.  More profit on each unit of A still sold.
2.  Fewer units of A sold.
3.  Some of the lost customers buy B, which the firm also owns.

At the profit-maximising price these effects balance. For brand A the
first-order condition is

$$s_A + (p_A - c_A)\frac{\partial s_A}{\partial p_A} + (p_B - c_B)\frac{\partial s_B}{\partial p_A} = 0.$$

The last term only appears because the firm owns both brands.

### Step 1: recover marginal costs

Costs are not in the data, but the first-order conditions reveal them.
Prices and shares are observed and the demand estimates give every
$\partial s / \partial p$, so the only unknowns are the costs. With one
equation per product, they can be solved for. The logic: if demand for a
brand is not very price sensitive and the firm still charges only a
little above cost, it would be leaving money on the table. Given the
elasticities, there is only one cost level at which the observed price
is the best one.

### Step 2: change who owns what

After the merger, the merged firm sets the prices of both firms’ brands.
Every brand of firm 1 now gets the third effect above for every brand of
firm 2, and vice versa. Raising the price of a firm 1 brand is less
costly than before, because customers who switch to a firm 2 brand are
no longer lost.

### Step 3: find the new equilibrium

Costs and demand stay the same. Only ownership changes. New prices solve
the first-order conditions with the new ownership. Rivals also respond:
when the merged firm raises prices, demand for rival brands rises, and
they raise their own prices a little.

### What the model leaves out

- **Efficiencies.** Mergers can lower costs. Here costs are held fixed,
  so the result is the price effect before any cost savings.
- **Coordination.** Firms are assumed to keep competing independently.
- **Entry, new products and repositioning** after the merger.

## 3. Data and demand

The PyBLP copy of the data includes firm identifiers. The price
coefficient is the IV logit estimate from part 1.

``` r
products <- read.csv("data/nevo_products.csv", stringsAsFactors = FALSE)
markets <- unique(products$market_ids)
market_rows <- split(seq_len(nrow(products)), factor(products$market_ids, levels = markets))
products$outside_share <- 1 - ave(products$shares, products$market_ids, FUN = sum)
products$delta <- log(products$shares) - log(products$outside_share)

products %>%
  group_by(Firm = firm_ids) %>%
  summarise(Brands = n_distinct(product_ids),
            `Average market share` = sum(shares) / length(markets)) %>%
  knitr::kable(digits = 3)
```

| Firm | Brands | Average market share |
|-----:|-------:|---------------------:|
|    1 |      9 |                0.199 |
|    2 |      9 |                0.158 |
|    3 |      2 |                0.049 |
|    4 |      3 |                0.051 |
|    6 |      1 |                0.020 |

Firms 1 and 2 are the two largest, with about 20% and 16% of the
potential market. They are the merger simulated here.

``` r
fit_linear <- function(y, X, Z = X) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  X_hat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(X_hat, X))
  setNames(as.vector(bread %*% crossprod(X_hat, y)), colnames(X))
}
iv_cols       <- grep("^demand_instruments", names(products), value = TRUE)
brand_dummies <- model.matrix(~ 0 + factor(product_ids), data = products)
alpha <- -unname(fit_linear(products$delta,
                            X = cbind(prices = products$prices, brand_dummies),
                            Z = cbind(as.matrix(products[iv_cols]), brand_dummies))["prices"])
alpha
```

    ## [1] 30.09776

## 4. The simulation

### Logit demand at any prices

In the logit, the derivatives and the new shares have formulas:

- $\partial s_j / \partial p_j = -\alpha s_j (1 - s_j)$ and
  $\partial s_j / \partial p_k = \alpha s_j s_k$.
- If prices change, mean utility changes by $-\alpha \Delta p$, and
  shares follow from the logit formula.

``` r
logit_shares <- function(delta) exp(delta) / (1 + sum(exp(delta)))

logit_derivatives <- function(s) {
  D <- alpha * outer(s, s)              # [j, k] = d s_j / d p_k
  diag(D) <- -alpha * s * (1 - s)
  D
}
```

### The ownership matrix

The first-order conditions of all products in a market can be written in
one line:

$$s + (\Omega \circ D')(p - c) = 0,$$

where $D'$ holds $\partial s_k / \partial p_j$ and $\Omega$ is the
ownership matrix: 1 if products $j$ and $k$ belong to the same firm, 0
otherwise. The $\circ$ multiplies element by element, so the zeros
switch off the third effect for brands owned by other firms. A merger is
a change in $\Omega$.

Solving for costs gives step 1:

$$c = p + (\Omega \circ D')^{-1} s.$$

``` r
ownership <- function(firms) outer(firms, firms, "==") * 1

marginal_costs <- function(p, s, D, own) p + as.vector(solve(own * t(D), s))
```

### Solving for post-merger prices

The first-order conditions can be rearranged as
$p = c - (\Omega \circ D(p)')^{-1} s(p)$: prices equal cost plus a
markup, and the markup depends on prices. The solution is found by
repetition. Start from the old prices, compute the markup, update
prices, and repeat until they stop changing.

Each update moves halfway towards the new prices. With full steps,
prices in two markets jump back and forth without settling; half steps
make every market converge.

``` r
simulate_merger <- function(p0, s0, firms, merging, damping = 0.5,
                            tol = 1e-12, max_iter = 1000) {
  own_pre  <- ownership(firms)
  firms_post <- ifelse(firms %in% merging, merging[1], firms)
  own_post <- ownership(firms_post)

  cost   <- marginal_costs(p0, s0, logit_derivatives(s0), own_pre)
  delta0 <- log(s0) - log(1 - sum(s0))

  p <- p0
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    s     <- logit_shares(delta0 - alpha * (p - p0))
    p_new <- cost - as.vector(solve(own_post * t(logit_derivatives(s)), s))
    if (max(abs(p_new - p)) < tol) {
      converged <- TRUE
      break
    }
    p <- p + damping * (p_new - p)
  }
  delta1 <- delta0 - alpha * (p - p0)

  list(price = p, share = logit_shares(delta1), cost = cost,
       converged = converged, iterations = iter,
       # Consumer surplus per person in the logit: log(1 + sum(exp(delta))) / alpha
       cs_change = (log(1 + sum(exp(delta1))) - log(1 + sum(exp(delta0)))) / alpha,
       spending  = sum(p0 * s0))
}
```

**Why does consumer surplus have that formula?** In the logit, the
expected utility of the best option has a closed form,
$\ln(1 + \sum_k e^{\delta_k})$. It is the same property of the
extreme-value distribution that gave the share formula. Dividing by the
price coefficient turns utility into money.

### Running it in every market

``` r
merging_firms <- c(1, 2)

results <- do.call(rbind, lapply(seq_along(markets), function(t) {
  rows <- market_rows[[t]]
  sim  <- simulate_merger(products$prices[rows], products$shares[rows],
                          products$firm_ids[rows], merging_firms)
  data.frame(market = markets[t],
             product = products$product_ids[rows],
             firm = products$firm_ids[rows],
             price_pre = products$prices[rows], price_post = sim$price,
             share_pre = products$shares[rows], share_post = sim$share,
             cost = sim$cost,
             converged = sim$converged, iterations = sim$iterations,
             cs_change = sim$cs_change, spending = sim$spending)
}))
results$merging      <- results$firm %in% merging_firms
results$price_change <- results$price_post / results$price_pre - 1
results$margin       <- (results$price_pre - results$cost) / results$price_pre
```

### Checks

``` r
c(markets_converged = sum(tapply(results$converged, results$market, all)),
  max_iterations = max(results$iterations),
  negative_costs = sum(results$cost < 0),
  median_margin = median(results$margin))
```

    ## markets_converged    max_iterations    negative_costs     median_margin 
    ##        94.0000000        37.0000000         1.0000000         0.3149889

**What these numbers say.**

- **Markets converged:** all 94 markets reached a stable set of prices,
  at most after 37 rounds.
- **Negative costs:** 1 out of 2256 products. A negative cost means the
  demand model cannot explain why the firm charges that price. One case
  in more than two thousand is negligible.
- **Median margin 31%:** before the merger, a typical brand’s price is
  about a third above its marginal cost. Part 1 showed that a firm
  selling a single product with an elasticity of about -3.7 would set a
  margin of about 27%. Margins here are higher because each firm owns
  several brands, and the customers who switch between its own brands
  are not lost.

## 5. Results

### Prices

``` r
results %>%
  group_by(Firm = firm) %>%
  summarise(Merging = ifelse(first(merging), "Yes", "No"),
            `Median price change (%)` = 100 * median(price_change),
            `Mean price change (%)` = 100 * mean(price_change),
            `Largest (%)` = 100 * max(price_change)) %>%
  knitr::kable(digits = 2)
```

| Firm | Merging | Median price change (%) | Mean price change (%) | Largest (%) |
|-----:|:--------|------------------------:|----------------------:|------------:|
|    1 | Yes     |                    5.40 |                  6.01 |       23.84 |
|    2 | Yes     |                    6.31 |                  7.52 |       40.84 |
|    3 | No      |                    0.07 |                  0.11 |        0.74 |
|    4 | No      |                    0.09 |                  0.13 |        0.64 |
|    6 | No      |                    0.02 |                  0.05 |        0.36 |

**What these numbers say.**

- **Merging firms:** the median brand of firm 1 rises 5.4% and of firm 2
  6.3%. This is the core result: the merger alone, with no change in
  costs or demand, makes these prices rise.
- **Largest increases:** in a few markets the increase is much larger,
  because the two firms’ brands there are large and the logit makes them
  close substitutes in proportion to their size. The median is the
  better summary.
- **Rivals:** prices rise by a fraction of a percent. When the merged
  firm raises prices, some of its customers turn to rivals, and rivals
  respond by raising their own prices a little.

### Quantities

``` r
quantities <- results %>%
  group_by(market) %>%
  summarise(merging_pre  = sum(share_pre[merging]),  merging_post = sum(share_post[merging]),
            rivals_pre   = sum(share_pre[!merging]), rivals_post  = sum(share_post[!merging]),
            outside_pre  = 1 - sum(share_pre),       outside_post = 1 - sum(share_post))

knitr::kable(data.frame(
  Group = c("Merging firms", "Rivals", "Outside good"),
  Change = c(sprintf("%.1f%% sales", 100 * mean(quantities$merging_post / quantities$merging_pre - 1)),
             sprintf("%.1f%% sales", 100 * mean(quantities$rivals_post / quantities$rivals_pre - 1)),
             sprintf("+%.1f percentage points of the market",
                     100 * mean(quantities$outside_post - quantities$outside_pre)))
))
```

| Group         | Change                               |
|:--------------|:-------------------------------------|
| Merging firms | -12.5% sales                         |
| Rivals        | 7.6% sales                           |
| Outside good  | +3.8 percentage points of the market |

### Consumers

``` r
welfare <- results %>%
  group_by(market) %>%
  summarise(cs_change = first(cs_change), spending = first(spending))

c(consumer_loss_as_pct_of_spending = -100 * mean(welfare$cs_change / welfare$spending))
```

    ## consumer_loss_as_pct_of_spending 
    ##                          4.21069

**What this number says.** For every 100 dollars consumers spent on
these cereals before the merger, they lose the equivalent of about 4.2
dollars.

The loss is expressed relative to what consumers spent on cereal before
the merger. It is somewhat smaller than the price increase times the
quantity bought before the merger, because consumers can avoid part of
the increase by switching to other brands or leaving the market.

## 6. Does the demand model matter?

Running the full simulation with the random-coefficients model requires
recomputing household-level shares at every candidate price, which is
beyond this project. Instead I compare the two models with the screening
measure competition authorities use first, the **Gross Upward Pricing
Pressure Index (GUPPI)**. It needs only pre-merger information:

$$\text{GUPPI}_j = \sum_{k \in \text{partner}} \text{Diversion}_{j \to k} \times \frac{p_k - c_k}{p_j}.$$

In words: of the customers brand $j$ loses when it raises its price, the
share that goes to the merging partner, times the margin earned on them,
relative to $j$’s price. It is the value of the sales that the merger
stops $j$ from losing. The higher it is, the stronger the incentive to
raise $j$’s price. Values of a few percent are usually treated as a sign
of possible concern.

Both inputs, diversion and margins, depend on the demand model.

``` r
guppi_market <- function(p, s, firms, D, merging) {
  cost <- marginal_costs(p, s, D, ownership(firms))
  div  <- -sweep(D, 2, diag(D), "/")   # [k, j]: share of j's lost customers that go to k
  diag(div) <- 0
  out <- data.frame(margin = (p - cost) / p, diversion_to_partner = NA_real_, guppi = NA_real_)
  for (j in which(firms %in% merging)) {
    partner <- which(firms %in% merging & firms != firms[j])
    out$diversion_to_partner[j] <- sum(div[partner, j])
    out$guppi[j] <- sum(div[partner, j] * (p[partner] - cost[partner])) / p[j]
  }
  out[firms %in% merging, ]
}
```

The random-coefficients derivatives come from part 2. That data uses
different product names, so each market is matched by position and the
prices are checked to be the same.

``` r
rc <- readRDS("output/rc_derivatives.rds")
stopifnot(length(rc) == length(markets))

guppi <- do.call(rbind, lapply(seq_along(markets), function(t) {
  rows  <- market_rows[[t]]
  p     <- products$prices[rows]
  s     <- products$shares[rows]
  firms <- products$firm_ids[rows]
  stopifnot(max(abs(rc[[t]]$price - p)) < 1e-6)
  rbind(cbind(Model = "Logit", guppi_market(p, s, firms, logit_derivatives(s), merging_firms)),
        cbind(Model = "Random coefficients", guppi_market(p, s, firms, rc[[t]]$derivs, merging_firms)))
}))

guppi_summary <- guppi %>%
  group_by(Model) %>%
  summarise(div = 100 * median(diversion_to_partner),
            margin = 100 * median(margin),
            guppi = 100 * median(guppi))

guppi_summary %>%
  rename(`Median diversion to partner (%)` = div,
         `Median margin (%)` = margin,
         `Median GUPPI (%)` = guppi) %>%
  knitr::kable(digits = 1)
```

| Model | Median diversion to partner (%) | Median margin (%) | Median GUPPI (%) |
|:---|---:|---:|---:|
| Logit | 17.0 | 32.9 | 5.6 |
| Random coefficients | 24.4 | 36.4 | 8.7 |

``` r
c(logit_simulated_median_price_increase = 100 * median(results$price_change[results$merging]))
```

    ## logit_simulated_median_price_increase 
    ##                              5.781919

**What these numbers say.**

- **Diversion to partner:** when a merging brand raises its price, the
  logit says about 17% of the customers it loses go to the other merging
  firm. Random coefficients say 24%, because the two firms own brands
  that are similar to each other.
- **Margin:** pre-merger margins are a little higher under random
  coefficients, because demand for each brand is estimated to be less
  sensitive to its own price.
- **GUPPI:** 5.6% under the logit against 8.7% under random
  coefficients. The logit value is close to the logit simulation’s
  median price increase of 5.8%, so the GUPPI is a reasonable shortcut
  for the size of the effect.

**Why this matters.** Under random coefficients, more of a brand’s lost
customers go to the merging partner, and margins are higher. Both raise
the GUPPI. This is the point of parts 1 and 2: the logit sends too many
customers out of the market and spreads the rest by market share, which
understates how closely the merging brands compete. With the richer
model the pricing pressure is larger, so the logit simulation is likely
a lower bound on the price effect.

## 7. Limitations

- The data are simulated and firm names are hidden, so the numbers
  illustrate the method and say nothing about a real merger.
- The full simulation uses the logit. The random-coefficients model is
  only compared through the GUPPI.
- Costs are fixed (no efficiencies) and firms are assumed to compete in
  prices before and after the merger.

## References

- Farrell, J., and Shapiro, C. (2010). Antitrust Evaluation of
  Horizontal Mergers: An Economic Alternative to Market Definition.
  *B.E. Journal of Theoretical Economics* 10(1).
- Nevo, A. (2000). Mergers with Differentiated Products: The Case of the
  Ready-to-Eat Cereal Industry. *RAND Journal of Economics* 31(3).
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy* 9(4).
- Werden, G., and Froeb, L. (1994). The Effects of Mergers in
  Differentiated Products Industries: Logit Demand and Merger Policy.
  *Journal of Law, Economics, & Organization* 10(2).
