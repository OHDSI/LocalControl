# Format and check code:
OhdsiRTools::formatRFolder()
OhdsiRTools::checkUsagePackage("LocalControl")
OhdsiRTools::updateCopyrightYearFolder()

# Create manual and vignettes:
shell("rm extras/LocalControl.pdf")
shell("R CMD Rd2pdf ./ --output=extras/LocalControl.pdf")
