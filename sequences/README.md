Sequence data
-------------

This directory is where the cluster-specific sequence alignments used by all downstream
analyses should go.  Due to GISAID restrictions, we are unable to include
the sequence files directly here.  Instead we provide information necessary to
necessary to construct these files yourself.

Firstly, the file `GISAID_Acknowledgement_Table.csv` lists the GISAID accession numbers
of all sequences used in the analyses. To access the corresponding sequences, you will
need to apply for an account at https://gisaid.org and to download them from there.
These sequences must then be aligned.

At this point, the individual sequences from the full alignment must be distributed between
the following fasta files:
- `australia.masked`
- `china.masked`
- `diamond_princess.masked`
- `dutch1.masked`
- `dutch2.masked`
- `french1.masked`
- `french2.masked`
- `iceland1.masked`
- `iceland2.masked`
- `iran.masked`
- `italy.masked`
- `spain.masked`
- `washington1.masked`
- `washington2.masked`
- `welsh.masked`

To assist with this the files *.masked.headers are provided, which are
just the header lines for each of the sequences which are to be added
to the corresponding *.masked file.  These header lines contain the
sequence accession numbers, and thus can be used to uniquely identify
which sequence belongs in each file.  When doing this, take care to
preserve the rest of the information in the header, as this is used in
downstream processing.
