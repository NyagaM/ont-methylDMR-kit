# Differential methylation analysis using bedmethyl files from long-read (ONT) data
`ont-methylDMR-kit` is a pipeline to call differentially methylated regions (DMRs) between haplotypes, two samples, or between two groups, coupled with annotation and visualisation, on long-read (ONT) sequencing bedmethyl files generated using [modkit](https://github.com/nanoporetech/modkit). This pipeline is inspired by my work on rare disorders, and the fact that long-read sequencing has the potential to comprehensively identify all modified bases, such as 5-Methylcytosine (5mC), 5-Hydroxymethylcytosine (5hmC), N6-methyladenine (6mA), and N4-methylcytosine (4mC) that have been identified for a growing number of rare disorders and imprinted disorders.

The pipeline is built using [Nextflow](https://www.nextflow.io/), a bioinformatics workflow manager that enables the development of portable and reproducible workflows. 
There are two docker images available from [DockerHub](https://hub.docker.com/repository/docker/nyagam/ont-methyl-kit/general) and [DockerHub](https://hub.docker.com/repository/docker/nyagam/modbamtools-v0.4.8/general) that contains all the tools/softwares required by the pipeline, making results highly reproducible. 

# DMR analysis:
DMR analysis is performed using the R-package [DSS](https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html). It supports calling DMRs across 5-Methylcytosine (`--5mC`), 5-Hydroxymethylcytosine (`--5hmC`), N6-methyladenine (`--6mA`), and N4-methylcytosine (`--4mC`) modified bases. It also supports haplotype-specific DMRs using `--phased_mC`, `--phased_mA`, and `--phased_hmC`. Just provide raw bedmethyl files generated using [modkit](https://github.com/nanoporetech/modkit) using either `--input_file1` and `--input_file2` for calling DMRs between haplotypes or two samples; or `--input_group1` and `--input_group2` for group analysis (bedmethyls in these two folders must have `*.bed` extension). The current default options for [DSS](https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html); delta (threshold for defining DMR) at 10%, p-values threshold for calling DMR at 0.01, minimum length (in basepairs) required for DMR methylation change analysis at 100, minimum number of CpG sites required for DMR at 10, and merging two DMRs that are very close to each other at less than 100 basepairs. See [DSS](https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html) for more information. 

# Annotation:
Significant DMRs are annotated to provide information on whether DMRs overlap with promoters, exons, and introns. A compressed file, annotations.zip (which needs to be unzipped `tar -xvf annotations.zip`), is provided with the pipeline that contains the annotation information, which is based on gencode v44. 

# Visualisation:
Annotated DMRs are plotted using [modbamtools](https://github.com/rrazaghi/modbamtools). It supports haplotype-specific DMR plotting (by providing haplotagged modified bam files using `--phased_modBam` ) or DMRs between two samples (by providing modified bam files using `--input_modBam1` and `--input_modBam2`) but does not support plotting DMRs from group analysis. You can provide a gene list with `--gene_list` to only plot DMRs for the provided genes.

# Installation and Usage:
```bash
$ git clone https://github.com/NyagaM/ont-methylDMR-kit.git
$ cd ont-methylDMR-kit
$ tar -xvf annotations.zip
```
To view usage and run options:

```bash
$ nextflow run main.nf --help
```
```bash
Usage: nextflow run main.nf [options]

Options:
  --input_file1     First input bedmethyl file from sample_1 or haplotype_1 bedmethyl (required)
  --input_file2     Second input bedmethyl file from sample_2 or haplotype_2 bedmethyl (required)
  --output_dir      Output directory (required)
  --input_modBam1   First input modified BAM used to generate bedmethyl for sample_1 (optional)
  --input_modBam2   Second input modified BAM used to generate bedmethyl for sample_2 (optional)
  --input_group1    Directory containing bedmethyl files (with *.bed extension) for group 1 (optional)
  --input_group2    Directory containing bedmethyl files (with *.bed extension) for group 2 (optional)
  --gene_list       A list of genes to generate plots on (optional)
  --5mC             Use this flag to trigger 5mC DMR calling (optional)
  --5hmC            Use this flag to trigger 5hmC DMR calling (optional)
  --6mA             Use this flag to trigger 6mA DMR calling (optional)
  --4mC             Use this flag to trigger 4mC DMR calling (optional)
  --phased_mC       Use this flag to trigger haplotagged 5mC/4mC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_mA       Use this flag to trigger haplotagged 6mA DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_hmC      Use this flag to trigger haplotagged 5hmC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_modBam   Haplotagged modified BAM (required for plotting DMRs if --phased_mC/--phased_mA is used)
  --help            Print this help message
```
