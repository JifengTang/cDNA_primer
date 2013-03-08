These scripts are for identifying potential full-length (FL) subreads or CCS reads using the 5' and 3' primer ligated to the transcripts during the cDNA library preparation step.

Refer to this wiki_ for detailed explanation on how to use the full-length identification scripts.

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



===========================================================                    
Extra filtering to eliminate subreads with missed adapters
===========================================================
If SMRTbell adapters are missed, sometimes it'll still be considered full-length by barcode_trimmer.py (especially
when the 5' and 3' primers are identical or highly similar). To further eliminate these subreads, after running
barcode_trimmer.py, you can run this extra script:

usage: filter_53seen.py TRIMMED_FASTA PRIMER_INFO OUTPUT_FILENAME

You should use the output fasta and .primer_info.txt from barcode_trimmer.py. 

Currently, this filtering is done by simply looking at the FL subread length distributions in the same ZMW
and eliminate those that have too short or too long subread length (despite seeing both 5' and 3').


The output from this script are all subreads that: (1) have 5' and 3' seen and (2) are likely not to contain a 
missing adapter.


===========================================================
Aligning subreads to known transcripts and plotting results
===========================================================
The alignQC.py is a lot more demanding on prerequisites. Read here_ to usage.

.. _here: https://github.com/Magdoll/cDNA_primer/wiki/Aligning-to-known-transcripts-for-QC


usage: alignQC.py [-h] -d OUTPUT_DIRECTORY -m PRIMER_MATCH_FILE -p OUTPUT_PREFIX [--read_pickle READ_PICKLE] [--ref_size REF_SIZE] [--refStrandPickle REFSTRANDPICKLE] [--restrictByPM] job_directory

Create some plots for transcript analyses.

positional arguments:
  job_directory

optional arguments:
  -h, --help            show this help message and exit
  -d OUTPUT_DIRECTORY   OUTPUT_DIRECTORY
  -m PRIMER_MATCH_FILE  PRIMER_MATCH_FILE
  -p OUTPUT_PREFIX      OUTPUT_PREFIX
  --read_pickle         READ_PICKLE
  --ref_size            REF_SIZE
  --refStrandPickle     REFSTRANDPICKLE
  --restrictByPM        Using .primer_info.txt to restrict what subreads to look at


