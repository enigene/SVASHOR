# SVASHOR
Set of scripts that generate reports with Structural Variants of Alpha Satellite
Higher-Order Repeats in the process of annotation human genome (GRCh38/hg38) assembly.

Before creating a report, you need to get a BED file annotated with HumAS-HMMER.

## Prerequisites
The pipeline requires the following components:
* awk, sed, and other standard unix command-line programs
* bedtools [bedtools utilities](http://bedtools.readthedocs.io)
* R [The R Project](https://www.r-project.org/) with these libraries:
    * xtable
    * xml2
    * stringr

## Simple usage
All BED files in the input directory get into processing:
```
svashor.sh /path/to/the/input/directory
```

## Output
Reports are obtained for all contigs in the region and printed as HTML files for
each HORs, also the directory `AS_HOR_structure` is created, which contains files
of an individual contigs and files with full HOR maps and statistics for the HOR
variants.

For example:
![HOR variants from chromosome 10](example/S1C10H1L.html)
