.PHONY: document install check

## regenerate man/ from the roxygen comments (requires roxygen2)
document:
	Rscript -e 'roxygen2::roxygenise()'

## install the package into the current R library
install:
	R CMD INSTALL .

## build the tarball and run a full R CMD check (this runs the test suite)
check:
	R CMD build . && R CMD check --no-manual "$$(ls -t saecol2026_*.tar.gz | head -n1)"
