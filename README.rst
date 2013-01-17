Some scripts for identifying potentially full-length transcripts
in PacBio transcriptome data.

This is work in progress. Scripts are not well tested :)
If you have problems please contact etseng@pacificbiosciences.com

Thanks!


Prerequisite for running the scripts
=====================================
1) Python 2.7.x
2) BioPython 

To identify 5'/3' primers you will need the PB Barcode module from PacBio's DevNet site:
http://www.smrtcommunity.com/Share/Code?id=a1q70000000GsyfAAC&strRecordTypeName=Code

download and follow directions. It will also ask that you install an old version of HMMER.

To get primer info for filtered_subreads.fasta::
    PacBioBarcodeIDCCS.py filtered_subreads.fasta primers.fa output

You will need to know the 5' and 3' primers used in the cDNA library prep.
The file format for primers.fa should be::
    >F0
    *5' sequence here*
    >R0
    *3' sequence here (but in reverse complement)*

If you have more than one set of primers you can add more, just name them F1,R1,F2,R2...

NOTE: if you have symmetrical 5'/3' primers, there's a chance my script for identifying polyA tail will fail. To prevent that one of the ways you can do is by adding a few As (like 4As or 5As) before the 3' primer sequence. In this way it's more likely that the 3' primer will be identified as the one following the polyA tail.

output/ is the output directory for PacBioBarcodeIDCCS.

In example/ I have put a test set that you can play with. 


Trimming away primers and polyA tails
========================================
To trim & output subreads that have BOTH 5' and 3' primers seen (if there is polyA will remove too)::
    scripts/barcode_trimmer.py -i <input_filename> -d <PBBarcode output directory> -o <output_filename>

For example::
    scripts/barcode_trimmer.py -i filtered_subreads.fasta -d output -o filtered_subreads.53seen_trimmed.fa


To output everything including no primers seen, use option --output-anyway.
To output all subreads but not necessarily 3', use --right-nosee-ok. Conversely, --left-no-seeok will output (and trim if necessary) all subreads that have 3' but not necessarily 5'.


See scripts/barcode_trimmer.py --help for the full set of parameters.

NOTE: ignore the primer percentages that barcode_trimmer.py outputs. It's not complete and to get accurate primer statistics use count_5seen.py below.


Explanation of .primer_info.txt
=================================
After you run barcode_trimmer.py you should get this .txt file which is a table where 5seen/polyAseen/3seen is '1' if it is seen, otherwise '0'. Sequence strand can be determined by observing 5' and/or 3' at either the beginning or end of the sequence.

You can use .primer_info.txt to get the accurate number of 5', 3', and 5'&3' primers seen per subreads and ZMW, do::
    scripts/count_5seen.py filtered_subreads.53seen_trimmed.fa.primer_info.txt filtered_subreads.fasta > filtered_subreads.53seen_trimmed.fa.primer_info.txt.summary

And filtered_subreads.53seen_trimmed.fa.primer_info.txt.summary will contain all the information.


Additional scripts that may be useful
=========================================
To get all filtered subreads that did not have CCS::
    scripts/grab_nonCCS_subreads.py <subreads_filename> <CCS_filename> <output_filename>

    For example::
        scripts/grab_nonCCS_subreads.py filtered_subreads.fasta filtered_CCS_subreads.fasta filtered_nonCCS_subreads.fasta

To sort a fasta file by length::
    scripts/sort_fasta_by_len.py <fasta_filename>

This will output a .sorted.fasta or .sorted.fa file.


