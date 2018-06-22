## Version 2.0.1 Release Notes
- Fixed bug with `KRp` and uneven sample size
- New warning issued when NAs present in the dataset
- Added Gamma distribution to `rsquared`
- Fixed bug to determine whether graph is cyclic
- Fixed bug thowing errors with calls to `lmer`
- Fixed bug creating errors with objects fit using `lmerTest`
- Added pkgdown website

## Version 2.0 Release Notes

### New syntax
-All functions have been re-written from the ground up
-Incorporates new `psem` function and S3 objects
-All necessary information can now be obtained with a single function `summary`

### Updated R[2] functions
-Extends to new distributions and model types using a single function `rsquared`

### New standardization procedures
-Implements range standardization for all response types
-Adds multiple forms of standardization for binary responses, see `?coefs`

### NOTES
-Removed `sem.plot` function
-Removed `sem.lavaan` function
-`groups=` argument is currently broken but will be fixed in version 2.1 (see doc)
