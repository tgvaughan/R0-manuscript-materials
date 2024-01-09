R_scripts
---------

This directory contains R scripts necessary for plotting the BEAST 2 results
and for performing auxilliary analyses.

Before running any of these scripts, make sure of the following:

1. The ../sequences directory must contain sequences downloaded from GISAID.
   For details on this, see the README.md in that directory.

2. The BEAST 2 XMLs in ../XMLs must be executed, leaving results in ../Results

3. The Bash script `extract_sequence_dates.sh` must be run.  (This extracts sequence
sample time information from the FASTA-formatted sequence files in
../sequences.)

Once these steps are complete, begin by running `estimateRe_clusters.R` to
perform an EpiEstim-based inference of R0 using sequence sample times only.

The full post-processing and figure generation script `make_plots.R` can then
be run.

Note that the `trajDataTools.R` is just a helper script required by `make_plots.R`
to parse the trajectory files produced by the BEAST 2 EpiInf package.

Generated figures will be left in the directory ../figures.
