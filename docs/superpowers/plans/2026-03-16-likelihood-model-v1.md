# likelihood.model v1.0.0 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Strip `likelihood.model` to v1.0.0, removing the contribution builder, named-distribution wrapper, and Weibull example. The package becomes the minimal, stable concept + inference foundation.

**Architecture:** Remove 3 R files (`model-contr.R`, `model-name.R`, `example-weibull.R`), drop R6 dependency, rewrite tests to use `exponential_lifetime` exclusively, update vignettes.

**Tech Stack:** R, S3 classes, roxygen2, testthat 3

**Spec:** `docs/superpowers/specs/2026-03-16-v1.0-and-likelihood.contr-design.md` (Package 1 section)

**Prerequisite:** `likelihood.contr` package must exist and pass tests first.

---

## File Structure

After v1.0.0, the R/ directory contains only:

```
R/
  core-generics.R        # concept, generics, defaults (kept, updated)
  core-fit.R             # fit, sampler (kept unchanged)
  core-fisher_mle.R      # fisher_mle, fisher_boot (kept unchanged)
  core-fisherian.R       # support, relative_likelihood, etc. (kept unchanged)
  core-lrt.R             # lrt (kept unchanged)
  example-exponential.R  # exponential_lifetime (kept unchanged)
```

Removed:
```
R/model-contr.R          # deleted
R/model-name.R           # deleted
R/example-weibull.R      # deleted
```

---

### Task 1: Remove source files and update DESCRIPTION

**Files:**
- Delete: `R/model-contr.R`
- Delete: `R/model-name.R`
- Delete: `R/example-weibull.R`
- Modify: `DESCRIPTION`
- Modify: `R/core-generics.R` (remove observed_info deprecation warning, keep as silent deprecated)

- [ ] **Step 1: Delete the three R source files**

```bash
git rm R/model-contr.R R/model-name.R R/example-weibull.R
```

- [ ] **Step 2: Update DESCRIPTION**

Change version to 1.0.0. Remove R6 from Imports. Update Description text to
remove mention of contribution models and named distributions. Keep it focused
on the concept and inference.

```
Version: 1.0.0
Imports:
    algebraic.mle (>= 2.0.0),
    generics,
    stats,
    numDeriv,
    boot
Suggests:
    mvtnorm,
    algebraic.dist,
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
```

Note: `dplyr` and `tibble` are dropped from Suggests (no longer used
after removing `likelihood_name` examples).

- [ ] **Step 3: Remove observed_info deprecation message**

In `R/core-generics.R`, keep `observed_info.likelihood_model` with the
`.Deprecated()` call. It stays exported for maskedcauses compatibility and
emits a deprecation warning so users know to migrate.

Update any tests that assert this warning: change from
`expect_warning(observed_info(model), "deprecated")` to wrapping
the call in `suppressWarnings()` if the test's purpose is to verify
the returned function (not the warning itself). Keep one explicit
deprecation test.

- [ ] **Step 4: Run devtools::document()**

```bash
Rscript -e 'devtools::document()' 2>&1
```

Verify that NAMESPACE no longer references R6, likelihood_contr_model,
likelihood_name, weibull_uncensored, likelihood_exact_weibull, or
prepare_args_list.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat!: strip to v1.0.0 - remove contr, name, weibull"
```

---

### Task 2: Rewrite tests

**Files:**
- Modify: `tests/test.R`

This is the largest task. Every test that references `likelihood_contr_model`,
`likelihood_name`, `weibull_uncensored`, or `likelihood_exact_weibull` must
be removed or rewritten to use `exponential_lifetime`.

- [ ] **Step 1: Identify tests to remove**

Search for all test blocks referencing removed features:

```bash
grep -n 'likelihood_contr_model\|likelihood_name\|weibull_uncensored\|likelihood_exact_weibull\|likelihood_name_model\|contr_model\|model-name\|observed_info' tests/test.R
```

Remove entire `test_that(...)` blocks that test removed features.

- [ ] **Step 2: Rewrite Fisherian inference tests**

Tests for `support`, `relative_likelihood`, `likelihood_interval`,
`profile_loglik`, and `evidence` that used `weibull_uncensored` must be
rewritten to use `exponential_lifetime`. Key pattern:

Before:
```r
model <- weibull_uncensored("x")
df <- data.frame(x = rweibull(50, shape = 2, scale = 1))
result <- fit(model)(df, par = c(shape = 1.5, scale = 0.8))
```

After:
```r
model <- exponential_lifetime("t")
df <- data.frame(t = rexp(50, rate = 2))
result <- fit(model)(df)  # closed-form, no initial guess needed
```

Note: `exponential_lifetime` is 1-parameter, so 2D profile likelihood tests
must be removed (there is no 2-parameter model in v1.0.0). The 1D versions
suffice.

- [ ] **Step 3: Rewrite LRT tests**

LRT tests that compare two different model types (e.g., exponential vs Weibull)
must be rewritten to compare two `exponential_lifetime` models with different
parameter values. This still tests the LRT machinery correctly.

- [ ] **Step 4: Remove FIM mock tests that are no longer needed**

The `mock_exp` FIM tests can stay (they test `fim.likelihood_model` with
a generic mock). Remove any that specifically test `weibull_uncensored` FIM.

- [ ] **Step 5: Run tests**

```bash
Rscript -e 'devtools::load_all(); testthat::test_file("tests/test.R")' 2>&1
```

Target: 0 failures. Test count will decrease (expected, since we removed
features).

- [ ] **Step 6: Commit**

```bash
git add tests/test.R
git commit -m "test: rewrite tests for v1.0.0 (exponential_lifetime only)"
```

---

### Task 3: Update vignettes

**Files:**
- Delete: `vignettes/likelihood-contributions.Rmd`
- Delete: `vignettes/likelihood-name-model.Rmd`
- Modify: `vignettes/getting-started.Rmd`
- Modify: `vignettes/algebraic-mle-integration.Rmd`

- [ ] **Step 1: Delete removed vignettes**

```bash
git rm vignettes/likelihood-contributions.Rmd
git rm vignettes/likelihood-name-model.Rmd
# Remove knitr artifacts if tracked; ignore errors for untracked files
git rm -f --ignore-unmatch vignettes/likelihood-contributions.R \
  vignettes/likelihood-contributions.html \
  vignettes/likelihood-name-model.R \
  vignettes/likelihood-name-model.html
```

- [ ] **Step 2: Rewrite getting-started.Rmd**

Rewrite to use only `exponential_lifetime`. Cover:
- Creating a model
- Fitting with `fit()`
- Examining results with `coef()`, `vcov()`, `confint()`, `summary()`
- Fisherian inference: `support()`, `likelihood_interval()`
- LRT between two models
- Bootstrap inference with `sampler()`

- [ ] **Step 3: Update algebraic-mle-integration.Rmd**

Replace any `weibull_uncensored` or `likelihood_name` examples with
`exponential_lifetime`. The vignette's algebraic.dist section (guarded
with `requireNamespace`) may need adjustment.

- [ ] **Step 4: Verify vignettes build**

```bash
Rscript -e 'devtools::build_vignettes()' 2>&1
```

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs: rewrite vignettes for v1.0.0"
```

---

### Task 4: Update NEWS.md and CLAUDE.md

**Files:**
- Modify: `NEWS.md`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add v1.0.0 section to NEWS.md**

Add a section documenting the breaking changes, removals, and migration path.
Reference `likelihood.contr` as the new home for contribution-based models.

- [ ] **Step 2: Update CLAUDE.md**

Remove sections about `likelihood_contr_model`, `likelihood_name`,
`weibull_uncensored`. Update the architecture section, file structure,
and dependency list. Add note about `likelihood.contr` as companion package.

- [ ] **Step 3: Commit**

```bash
git add NEWS.md CLAUDE.md
git commit -m "docs: update NEWS.md and CLAUDE.md for v1.0.0"
```

---

### Task 5: Full verification

- [ ] **Step 1: Run devtools::document()**

```bash
Rscript -e 'devtools::document()' 2>&1
```

- [ ] **Step 2: Run R CMD check**

```bash
Rscript -e 'devtools::check(args = c("--as-cran", "--no-manual"))' 2>&1
```

Target: 0 errors, 0 warnings.

- [ ] **Step 3: Run coverage**

```bash
Rscript -e 'cov <- covr::package_coverage(); cat(sprintf("Coverage: %.1f%%\n", covr::percent_coverage(cov)))' 2>&1
```

Target: >= 95%.

- [ ] **Step 4: Verify downstream packages**

```bash
# Each should install and pass tests against the stripped likelihood.model
Rscript -e 'devtools::install(); devtools::check("/home/spinoza/github/rlang/maskedcauses")' 2>&1
Rscript -e 'devtools::check("/home/spinoza/github/rlang/flexhaz")' 2>&1
Rscript -e 'devtools::check("/home/spinoza/github/rlang/maskedhaz")' 2>&1
```

Target: 0 errors in all three.

- [ ] **Step 5: Verify likelihood.contr against stripped v1.0.0**

```bash
Rscript -e 'devtools::install(); devtools::check("/home/spinoza/github/rlang/likelihood.contr")' 2>&1
```

Target: 0 errors.

- [ ] **Step 6: Commit any fixes**

```bash
git add -A
git commit -m "chore: pass full verification pipeline for v1.0.0"
```
