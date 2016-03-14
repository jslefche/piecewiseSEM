## Update package
This is a feature and bug-fix update that adds new plotting capabilities, fixes long-standing issues with `lme4` models, and issues with construction of goodness-of-fit tests.

## R CMD check results
There were no ERRORs or WARNINGs. 

## NOTES 
* I have added "par" and "legend" to the namespace based on returned comments

* Several examples take longer >5 to run on the test system but are crucial for reproducing the analyses in the primary literature. I have commented out several of the examples in an effort to reduce the run time

## Test environments
* local OS X install, R 3.2.3
* win-builder (R-devel and R-release)