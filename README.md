LocalControl
============

Introduction
============

LocalControl is an R package implementing novel nonparametric approaches to address biases and confounding when comparing treatments or exposures in observational studies of outcomes.

See R help pages for examples on how to run each of the package functions.

Getting Involved
================
* Package manual: [LocalControl manual](https://CRAN.R-project.org/package=LocalControl/LocalControl.pdf)
* Package vignette: [LocalControl vignette](https://CRAN.R-project.org/package=LocalControl/vignettes/LocalControl-jss-2020.pdf)
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="https://github.com/OHDSI/LocalControl/issues">GitHub issue tracker</a> for all bugs/issues/enhancements

Building considerations for CRAN
================================
Given we have a large manual, building a compressed pdf of the vignette may be required to not have the submission rejected.
Within RStudio this can be done via:

```r
devtools::build(args = c('--compact-vignettes=both'))
```

License
=======
LocalControl is licensed under Apache License 2.0
