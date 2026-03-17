## R CMD check results

0 errors | 0 warnings | 0 notes

## Update (v0.9.1 → v1.0.0)

Removed re-export chain for algebraic.mle generics. Now uses Depends:
algebraic.mle and imports algebraic.dist directly for S3 method registration.

## Coordinated submission

This is part of a coordinated 6-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- algebraic.dist 1.0.0
- algebraic.mle 2.0.2
- likelihood.model 1.0.0 (this package)
- compositional.mle 2.0.0
- flexhaz 0.5.1
- maskedcauses 0.9.3

## Test environments

* local Ubuntu 24.04, R 4.3.3
