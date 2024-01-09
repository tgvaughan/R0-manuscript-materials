BEAST 2 XMLs
------------

This directory contains the BEAST 2 XML files used to produce the
Bayesian phylodynamic estimates of R0 and total case count.

Before running these scripts, firstly ensure the GISAID sequence data
has been downloaded and is provided in the appropriate format in the
../sequences directory, as described in the README.md file found in
that directory.

Secondly, ensure you have BEAST 2 installed (and accessible from the
command line), and that the following BEAST 2 packages are installed:
* BDSKY (Used for birth-death skyline plot inference)
* EpiInf (Used for prevalence trajectory incidence)
* feast (General-purpose library)

To actually run the analyses, execute the script `run_analyses.sh` from
the command line in this directory.

**Warning** These analyses can take a long time to finish.  We
recommend you run them on a dedicated cluster.
