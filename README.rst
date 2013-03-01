These scripts are for identifying potential full-length (FL) subreads or CCS reads using the 5' and 3' primer ligated to the transcripts during the cDNA library preparation step.

Refer to this wiki_ for detailed explanation on how to use the scripts in this repository.

.. _wiki: https://github.com/Magdoll/cDNA_primer/wiki/How-to-identify-full-length-transcripts-in-PacBio-data



usage: barcode_trimmer.py [-h] -i INPUT_FASTA -d BARCODE_REPORT_DIR -o
                          OUTPUT_FILENAME [--left-nosee-ok] [--right-nosee-ok]
                          [--output-anyway] [--change-seqid]
                          [--min-seqlen MIN_SEQLEN]

Trim barcode

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA        Input fasta filename (default: None)
  -d BARCODE_REPORT_DIR       Barcode report directory (default: None)
  -o OUTPUT_FILENAME    Output filename (default: None)
  --left-nosee-ok       OK if 5' end not detected (default: False)
  --right-nosee-ok      OK if 3' end not detected (default: False)
  --output-anyway       Still output seqs w/ no barcode (default: False)
  --change-seqid        Change subread id to reflect trimming (default: False)
  --min-seqlen MIN_SEQLEN
                        Minimum seqlength to output (default 50) (default: 50)



Extra filtering to eliminate subreads with missed adapters.
