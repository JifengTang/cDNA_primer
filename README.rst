Some scripts for identifying potentially full-length transcripts
in PacBio transcriptome data.

This is work in progress. Scripts are not well tested :)
If you have problems please contact etseng@pacificbiosciences.com

Thanks!

File explanation
===================
XXX.53seen_trimmed.fa --- filtered subhead sequences that have both 5' and 3' primer seen, with the primers (and polyA tail) trimmed.
XXX.53seen_trimmed.primer_info.txt --- a table where 5seen/polyAseen/3seen is '1' if it is seen, otherwise '0'. Sequence strand can be determined by observing 5' and/or 3' at either the beginning or end of the sequence.

If you want to get sequences that just have 5' or 3' seen, see below.


Prerequisite for running the scripts
=====================================
1) Python 2.7.x
2) BioPython 

To identify 5'/3' primers you will need the PB Barcode module from PacBio's DevNet site:
http://www.smrtcommunity.com/Share/Code?id=a1q70000000GsyfAAC&strRecordTypeName=Code

download and follow directions. It will also ask that you install an old version of HMMER.

For example I run the following command to get primer info for CCS reads::
    PacBioBarcodeIDCCS.py filtered_CCS_subreads.fasta primers.fa outputCCS
    PacBioBarcodeIDCCS.py filtered_subreads.fasta primers.fa output

You will need to know the 5' and 3' primers used in the cDNA library prep.
The file format for primers.fa should be:

>F0
5' sequence here
>R0
3' sequence here (but in reverse complement)

And if you have more than one set of primers you can add more:
>F1
>R1
....


output/ is the output directory for PacBioBarcodeIDCCS.

In example/ I have put a test set that you can play with. 


Commands
==================
To get all filtered subreads that did not have CCS::
    scripts/grab_nonCCS_subreads.py <subreads_filename> <CCS_filename> <output_filename>

For example::
    scripts/grab_nonCCS_subreads.py filtered_subreads.fasta filtered_CCS_subreads.fasta filtered_nonCCS_subreads.fasta

  

To trim & output subreads that have BOTH 5' and 3' primers seen (if there is polyA will remove too)::
    scripts/barcode_trimmer.py -i <input_filename> -d <PBBarcode output directory> -o <output_filename>

For example::
    scripts/barcode_trimmer.py -i filtered_nonCCS_subreads.fasta -d output -o filtered_nonCCS_subreads.53seen_trimmed.fa 
    scripts/barcode_trimmer.py -i filtered_CCS_subreads.fasta -d outputCCS -o filtered_CCS_subreads.53seen_trimmed.fa


To output everything including no primers seen, use option --output-anyway.
To output all subreads but not necessarily 3', use --right-nosee-ok. Conversely, --left-no-seeok will output (and trim if necessary) all subreads that have 3' but not necessarily 5'.


See scripts/barcode_trimmer.py --help for the full set of parameters.

NOTE: ignore the primer percentages that barcode_trimmer.py outputs. It's not complete and to get accurate primer statistics use count_5seen.py below.


Finally, to get the accurate number of 5', 3', and 5'&3' primers seen per subreads and ZMW, do::
    scripts/count_5seen.py filtered_CCS_subreads.53seen_trimmed.fa.primer_info.txt filtered_subreads.fasta > filtered_CCS_subreads.53seen_trimmed.fa.primer_info.txt.summary

And filtered_CCS_subreads.53seen_trimmed.fa.primer_info.txt.summary will contain all the information.


To sort a fasta file by length::
    scripts/sort_fasta_by_len.py <fasta_filename>

This will output a .sorted.fasta or .sorted.fa file.


