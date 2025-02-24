#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Help message
def helpMessage = """
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
  --plots_only      Only run the plotting processes, requires --annotated_dmrs and BAM files
  --annotated_dmrs  Path to annotated DMR bed file (required if --plots_only is used)
  --help            Print this help message
"""

params.help = false
params.input_modBam1 = false
params.input_modBam2 = false
params.gene_list = false
params.phased_modBam = ""
params.input_file1 = ""
params.input_file2 = ""
params.input_group1 = ""
params.input_group2 = ""
params.output_dir = ""
params.'5mC' = false
params.'5hmC' = false
params.'6mA' = false
params.'4mC' = false
params.'phased_mA' = false
params.'phased_mC' = false
params.'phased_hmC' = false
params.plots_only = false
params.annotated_dmrs = ""

// Help information
if (params.help) {
    println helpMessage
    exit 0
}

// Separate validation for plots_only mode
if (params.plots_only) {
    if (!params.annotated_dmrs || !params.output_dir) {
        println "In plots-only mode, --annotated_dmrs and --output_dir are required"
        println helpMessage
        exit 0
    }
} else {
    // Validation for regular analysis mode
    if ((!params.input_file1 && !params.input_group1) || 
        (!params.input_file2 && !params.input_group2) || 
        !params.output_dir) {
        println "In analysis mode, input files (--input_file1/2 or --input_group1/2) and --output_dir are required"
        println helpMessage
        exit 0
    }
}

// Check if output directory exists
def outputDir = file(params.output_dir)
if (outputDir.exists()) {
    println "Output directory exists: ${params.output_dir}"
} else {
    println "Output directory does not exist. Creating one: ${params.output_dir}"
    outputDir.mkdirs()
}

// Input channels
input_ch1 = params.input_file1 ? Channel.fromPath(params.input_file1, checkIfExists: true) : Channel.empty()
input_ch2 = params.input_file2 ? Channel.fromPath(params.input_file2, checkIfExists: true) : Channel.empty()
annotationFile = file("${workflow.projectDir}/annotations/gencode.v44.annotation.exon-promoters-introns.sorted.bed")
modbamtools_gtf = file("${workflow.projectDir}/annotations/modbamtools.annotation.sorted.gtf.gz")
annotated_dmrs_ch = params.plots_only ? Channel.fromPath(params.annotated_dmrs, checkIfExists: true) : Channel.empty()

// Check for BAM files
bam_files_provided = params.input_modBam1 && params.input_modBam2
phased_bam_provided = (params.'phased_mC' || params.'phased_mA') && params.phased_modBam

input_bam_ch1 = bam_files_provided ? Channel.fromPath(params.input_modBam1, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()
input_bam_ch2 = bam_files_provided ? Channel.fromPath(params.input_modBam2, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()
phased_bam_ch = phased_bam_provided ? Channel.fromPath(params.phased_modBam, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()

// Create an empty gene list file if none is provided
empty_gene_list = file("${params.output_dir}/empty_gene_list.txt")
if (!params.gene_list) {
    empty_gene_list.text = ""
}

// Process 5mC methylbed file sample 1
process prep_bedmethyl_5mC_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 2' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 5mC methylbed file sample 2
process prep_bedmethyl_5mC_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 2' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 5-hydroxymethylcytosine (5hmC) methylbed file sample 1
process prep_bedmethyl_5hmC_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "h"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 5-hydroxymethylcytosine (5hmC) methylbed file sample 2
process prep_bedmethyl_5hmC_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "h"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 6mA methylbed file sample 1
process prep_bedmethyl_6mA_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "a"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 6mA methylbed file sample 2
process prep_bedmethyl_6mA_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "a"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 4mC methylbed file sample 1
process prep_bedmethyl_4mC_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Process 4mC methylbed file sample 2
process prep_bedmethyl_4mC_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged 5-hydroxymethylcytosine (5hmC) file sample 1
process prep_bedmethyl_hmC_HP_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "h"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged  5-hydroxymethylcytosine (5hmC) file sample 2
process prep_bedmethyl_hmC_HP_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "h"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged methylcytosine (5mC) file sample 1
process prep_bedmethyl_mC_HP_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged methylcytosine (5mC) file sample 2
process prep_bedmethyl_mC_HP_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "m"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged 6mA methylbed file sample 1
process prep_bedmethyl_6mA_HP_1 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "a"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Processing haplotagged 6mA methylbed file sample 2
process prep_bedmethyl_6mA_HP_2 {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_modified.bed"), emit: modified_bed
  script:
  """
  sed -i.bak 's/ /\\t/g' ${bed}
  awk -F'\\t' '\$4 == "a"' ${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > ${bed.baseName}_modified.bed
  """
}

// Evaluate significant DMRs using DSS
process dmr_calling {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/DMRs", mode: 'copy'
  cpus 12
  memory '16 GB'
  time 2.hours

  input:
    path bed1
    path bed2

  output:
    path("dmrs_table.bed"), emit: dmr_beds
    path("dmr_status.log"), emit: status_log
  script:
  """
  # Header for bed files
  echo -e "chr\\tpos\\tN\\tX" > header.txt
  # Initialize an empty file for merged results
  echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed

  # Determine available chromosomes in the input BED files
  available_chromosomes=\$(awk '{print \$1}' ${bed1} ${bed2} | sort -k1,1V | uniq | grep -E '^(chr)?[0-9X]+\$')

  if [ -z "\$available_chromosomes" ]; then
    echo "No valid chromosomes found in the input BED files." > dmr_status.log
    exit 0
  fi

  echo -e "Found and processing chromosomes:\n\$available_chromosomes" > dmr_status.log
  
  # Process each available chromosome
  for chr in \$available_chromosomes; do
    echo "Processing chromosome: \${chr}" >&2

    awk -v chr="\${chr}" '\$1 == chr' ${bed1} > ${bed1.baseName}_\${chr}.tmp.bed
    cat header.txt ${bed1.baseName}_\${chr}.tmp.bed > ${bed1.baseName}_\${chr}.bed
    rm ${bed1.baseName}_\${chr}.tmp.bed

    awk -v chr="\${chr}" '\$1 == chr' ${bed2} > ${bed2.baseName}_\${chr}.tmp.bed
    cat header.txt ${bed2.baseName}_\${chr}.tmp.bed > ${bed2.baseName}_\${chr}.bed
    rm ${bed2.baseName}_\${chr}.tmp.bed

    cat <<EOF > ont-methyl-kit_\${chr}.R
    library(DSS)
    require(bsseq)
    dat1 = read.table("${bed1.baseName}_\${chr}.bed", header=TRUE)
    dat2 = read.table("${bed2.baseName}_\${chr}.bed", header=TRUE)
    BSobj = makeBSseqData( list(dat1, dat2), c("C1","C2") )
    dmlTest = DMLtest(BSobj, group1=("C1"), group2=("C2"), smoothing=TRUE, ncores=12)
    dmrs = callDMR(dmlTest, delta=0.10, p.threshold=0.01, minlen=100, minCG=10, dis.merge=100, pct.sig=0.5)
    write.table(dmrs, file="dmrs_table_\${chr}.txt", row.names=FALSE, quote=FALSE, sep="\\t")
EOF

    Rscript ont-methyl-kit_\${chr}.R

    # Append the results to the merged file
    if [[ -f dmrs_table_\${chr}.txt ]]; then
      tail -n +2 dmrs_table_\${chr}.txt >> dmrs_table.bed
    else
      echo "Warning: dmrs_table_\${chr}.txt not found, skipping this chromosome"
    fi
  done
  """
}

// Evaluate significant group DMRs using DSS
process group_dmr_calling_5mC {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/group_DMRs", mode: 'copy'
  cpus 12
  memory '36 GB'
  time 6.hours

  input:
    path group1_beds
    path group2_beds

  output:
    path "dmrs_table.bed", emit: dmr_beds
    path "debug", emit: debug_files

  script:
  """
    # Create directories for group1 and group2 and debug
    mkdir -p group1 group2 debug

    # Copy and prep bed files to respective directories
    for bed in ${group1_beds}; do
      base_name=\$(basename \${bed} .bed)
      cp \${bed} group1/
      mv group1/\${bed} group1/\${base_name}.group1.bed
    done

    for bed in ${group2_beds}; do
      base_name=\$(basename \${bed} .bed)
      cp \${bed} group2/
      mv group2/\${bed} group2/\${base_name}.group2.bed
    done

    # Create header files
    echo -e "chr\\tpos\\tN\\tX" > header.tsv
    echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed

    # Process group1 bed files
    for bed in group1/*.bed; do
      available_chromosomes=\$(awk '{print \$1}' \${bed} | sort -k1,1V | uniq | grep -E '^(chr)?[0-9X]+\$')
      echo -e "\$available_chromosomes" > dmr_status.log
      for chr in \$available_chromosomes; do
        echo "Processing chromosome: \${chr}" >&2
        mkdir -p group1/\${chr}
        base_name=\$(basename \${bed} .group1.bed)
        sed -i.bak 's/ /\\t/g' \${bed}
        awk -F'\\t' '\$4 == "m"' \${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > group1/\${base_name}.group1.modified.bed
        awk -v chr="\${chr}" '\$1 == chr' group1/\${base_name}.group1.modified.bed > group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed
        cat header.tsv group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed > group1/\${chr}/\${base_name}.group1.modified_\${chr}.bed
        rm group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed
      done  
    done

    # Process group2 bed files
    for bed in group2/*.bed; do
      available_chromosomes=\$(awk '{print \$1}' \${bed} | sort -k1,1V | uniq | grep -E '^(chr)?[0-9X]+\$')
      echo -e "Found and processing chromosomes:\n\$available_chromosomes" > dmr_status.log
      for chr in \$available_chromosomes; do
        echo "Processing chromosome: \${chr}" >&2
        mkdir -p group2/\${chr}
        base_name=\$(basename \${bed} .group2.bed)
        sed -i.bak 's/ /\\t/g' \${bed}
        awk -F'\\t' '\$4 == "m"' \${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > group2/\${base_name}.group2.modified.bed
        awk -v chr="\${chr}" '\$1 == chr' group2/\${base_name}.group2.modified.bed > group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed
        cat header.tsv group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed > group2/\${chr}/\${base_name}.group2.modified_\${chr}.bed
        rm group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed
      done  
    done

    # Run DMR analysis per chrom
    available_chromosomes=\$(awk '{print \$1}' dmr_status.log | grep -E '^(chr)?[0-9X]+\$') && echo -e \${available_chromosomes}
    for chr in \$available_chromosomes; do
      cat <<EOF > ont-methyl-kit_grouped_\${chr}.R
      library(DSS)
      require(bsseq)
      dat1_files = Sys.glob("group1/\${chr}/*_\${chr}.bed")
      dat2_files = Sys.glob("group2/\${chr}/*_\${chr}.bed")
      if (length(dat1_files) == 0 || length(dat2_files) == 0) {
        stop("No files found for chromosome \${chr} in one of the groups.")
      }
      writeLines(dat1_files, con="dat1_files_\${chr}.txt")
      writeLines(dat2_files, con="dat2_files_\${chr}.txt")
      dat1_list = lapply(dat1_files, read.table, header=TRUE)
      dat2_list = lapply(dat2_files, read.table, header=TRUE)
      if (length(dat1_list) == 0 || length(dat2_list) == 0) {
        stop("No data read for chromosome \${chr} in one of the groups.")
      }
      BSobj = makeBSseqData(c(dat1_list, dat2_list), c(paste0("C", 1:length(dat1_list)), paste0("N", 1:length(dat2_list))))
      dmlTest = DMLtest(BSobj, group1=paste0("C", 1:length(dat1_list)), group2=paste0("N", 1:length(dat2_list)), smoothing=TRUE, ncores=12)
      dmrs = callDMR(dmlTest, delta=0.10, p.threshold=0.01, minlen=100, minCG=10, dis.merge=100, pct.sig=0.5)
      write.table(dmrs, file="dmrs_table_\${chr}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
EOF
      Rscript ont-methyl-kit_grouped_\${chr}.R

      # Append the results to the merged file
      if [[ -f dmrs_table_\${chr}.tsv ]]; then
        tail -n +2 dmrs_table_\${chr}.tsv >> dmrs_table.bed
      else
        echo "Warning: dmrs_table_\${chr}.tsv not found, skipping this chromosome"
      fi
    done
    mv *.txt *.R debug && rm -r group1 group2
  """
}

// Evaluate significant group DMRs using DSS
process group_dmr_calling_6mA {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/group_DMRs", mode: 'copy'
  cpus 12
  memory '36 GB'
  time 4.hours

  input:
    path group1_beds
    path group2_beds

  output:
    path "dmrs_table.bed", emit: dmr_beds
    path "debug", emit: debug_files

  script:
  """
    # Create directories for group1 and group2 and debug
    mkdir -p group1 group2 debug

    # Copy and prep bed files to respective directories
    for bed in ${group1_beds}; do
      base_name=\$(basename \${bed} .bed)
      cp \${bed} group1/
      mv group1/\${bed} group1/\${base_name}.group1.bed
    done

    for bed in ${group2_beds}; do
      base_name=\$(basename \${bed} .bed)
      cp \${bed} group2/
      mv group2/\${bed} group2/\${base_name}.group2.bed
    done

    # Create header files
    echo -e "chr\\tpos\\tN\\tX" > header.tsv
    echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed

    # Process group1 bed files
    for bed in group1/*.bed; do
      available_chromosomes=\$(awk '{print \$1}' \${bed} | sort -k1,1V | uniq | grep -E '^(chr)?[0-9X]+\$')
      echo -e "\$available_chromosomes" > dmr_status.log
      for chr in \$available_chromosomes; do
        echo "Processing chromosome: \${chr}" >&2
        mkdir -p group1/\${chr}
        base_name=\$(basename \${bed} .group1.bed)
        sed -i.bak 's/ /\\t/g' \${bed}
        awk -F'\\t' '\$4 == "a"' \${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > group1/\${base_name}.group1.modified.bed
        awk -v chr="\${chr}" '\$1 == chr' group1/\${base_name}.group1.modified.bed > group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed
        cat header.tsv group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed > group1/\${chr}/\${base_name}.group1.modified_\${chr}.bed
        rm group1/\${chr}/\${base_name}.group1.modified_\${chr}.tmp.bed
      done  
    done

    # Process group2 bed files
    for bed in group2/*.bed; do
      available_chromosomes=\$(awk '{print \$1}' \${bed} | sort -k1,1V | uniq | grep -E '^(chr)?[0-9X]+\$')
      echo -e "Found and processing chromosomes:\n\$available_chromosomes" > dmr_status.log
      for chr in \$available_chromosomes; do
        echo "Processing chromosome: \${chr}" >&2
        mkdir -p group2/\${chr}
        base_name=\$(basename \${bed} .group2.bed)
        sed -i.bak 's/ /\\t/g' \${bed}
        awk -F'\\t' '\$4 == "a"' \${bed} | awk -F'\\t' '\$5 >= 5' | awk -F'\\t' '{print \$1"\\t"\$3"\\t"\$5"\\t"\$12}' > group2/\${base_name}.group2.modified.bed
        awk -v chr="\${chr}" '\$1 == chr' group2/\${base_name}.group2.modified.bed > group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed
        cat header.tsv group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed > group2/\${chr}/\${base_name}.group2.modified_\${chr}.bed
        rm group2/\${chr}/\${base_name}.group2.modified_\${chr}.tmp.bed
      done  
    done

    # Run DMR analysis per chrom
    available_chromosomes=\$(awk '{print \$1}' dmr_status.log | grep -E '^(chr)?[0-9X]+\$') && echo -e \${available_chromosomes}
    for chr in \$available_chromosomes; do
      cat <<EOF > ont-methyl-kit_grouped_\${chr}.R
      library(DSS)
      require(bsseq)
      dat1_files = Sys.glob("group1/\${chr}/*_\${chr}.bed")
      dat2_files = Sys.glob("group2/\${chr}/*_\${chr}.bed")
      if (length(dat1_files) == 0 || length(dat2_files) == 0) {
        stop("No files found for chromosome \${chr} in one of the groups.")
      }
      writeLines(dat1_files, con="dat1_files_\${chr}.txt")
      writeLines(dat2_files, con="dat2_files_\${chr}.txt")
      dat1_list = lapply(dat1_files, read.table, header=TRUE)
      dat2_list = lapply(dat2_files, read.table, header=TRUE)
      if (length(dat1_list) == 0 || length(dat2_list) == 0) {
        stop("No data read for chromosome \${chr} in one of the groups.")
      }
      BSobj = makeBSseqData(c(dat1_list, dat2_list), c(paste0("C", 1:length(dat1_list)), paste0("N", 1:length(dat2_list))))
      dmlTest = DMLtest(BSobj, group1=paste0("C", 1:length(dat1_list)), group2=paste0("N", 1:length(dat2_list)), smoothing=TRUE, ncores=12)
      dmrs = callDMR(dmlTest, delta=0.10, p.threshold=0.01, minlen=100, minCG=10, dis.merge=100, pct.sig=0.5)
      write.table(dmrs, file="dmrs_table_\${chr}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
EOF
      Rscript ont-methyl-kit_grouped_\${chr}.R

      # Append the results to the merged file
      if [[ -f dmrs_table_\${chr}.tsv ]]; then
        tail -n +2 dmrs_table_\${chr}.tsv >> dmrs_table.bed
      else
        echo "Warning: dmrs_table_\${chr}.tsv not found, skipping this chromosome"
      fi
    done
    mv *.txt *.R debug && rm -r group1 group2

  """
}

// Annotate significant DMRs
process annotate_dmrs {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/annotated_DMRs", mode: 'copy'
  cpus 2
  memory '2 GB'

  input:
    path bed
    path annotationFile

  output:
    path("dmrs_table_annotated.bed"), emit: annotated_dmr_beds
  script:
  """
  echo 'chr\tstart\tend\tlength\tnSites\tmeanMethy1\tmeanMethy2\tdiff.Methy\tareaStat\tannotation_chr\tannotation_start\tannotation_end\tstrand\tannotation\tbiotype\tgene' > annotation_header.txt
  bedtools intersect -a dmrs_table.bed -b ${annotationFile} -wa -wb > dmrs_table_annotated.tmp
  cat annotation_header.txt dmrs_table_annotated.tmp > dmrs_table_annotated.bed
  """
}

// Plot significant DMRs using modbamtools
process plot_dmrs {
  label 'plot_methyl_dmrs'
  publishDir "${params.output_dir}/annotated_DMRs/plots", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    path gtf_file
    path gene_list
    tuple val(name1), file(bam1), file(bai1)
    tuple val(name2), file(bam2), file(bai2)

  output:
    path("*.html"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs
  script:
  """
  gene_list=${params.gene_list ? params.gene_list : empty_gene_list}
  genes_of_interest="genes_of_interest.txt"
  cp \${gene_list} \${genes_of_interest}

  # Create an empty error log file
  touch error.log

  plot_count=0

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nSites meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s genes_of_interest.txt ]] || grep -q -w "\${gene}" \${genes_of_interest}; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      modbamtools plot -r \${region} -g ${modbamtools_gtf} -s ${bam1.baseName},${bam2.baseName} -p \${output_prefix} ${bam1} ${bam2} -o ./ || {
        echo "Error processing region \${region} for gene \${gene}" >> error.log
      }
      plot_count=\$((plot_count + 1))
    fi
  done

  if [[ \$plot_count -eq 0 ]]; then
    echo "No plots generated, exiting successfully."
    exit 0
  fi
  """
}

// Plot haplotagged DMRs using modbamtools
process plot_haplotagged_dmrs {
  label 'plot_methyl_dmrs'
  publishDir "${params.output_dir}/annotated_DMRs/haplotagged_DMR_plots", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    path gtf_file
    path gene_list
    tuple val(name), file(bam), file(bai)
    
  output:
    path("*.html"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs
  script:
  """
  gene_list=${params.gene_list ? params.gene_list : empty_gene_list}
  genes_of_interest="genes_of_interest.txt"
  cp \${gene_list} \${genes_of_interest}

  # Create an empty error log file
  touch error.log
  plot_count=0

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s genes_of_interest.txt ]] || grep -q -w "\${gene}" \${genes_of_interest}; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      modbamtools plot -r \${region} -g ${modbamtools_gtf} -s ${bam.baseName} -p \${output_prefix} -hp ${bam} -o ./ || {
        echo "Error processing region \${region} for gene \${gene}" >> error.log
      }
      plot_count=\$((plot_count + 1))
    fi
  done

  if [[ \$plot_count -eq 0 ]]; then
    echo "No plots generated, exiting successfully."
    exit 0
  fi
  """
}

workflow {
    if (params.plots_only) {
        println "Running in plots-only mode..."
        if (params.'phased_mC' || params.'phased_mA') {
            if (phased_bam_provided) {
                plot_haplotagged_dmrs(annotated_dmrs_ch, 
                                    modbamtools_gtf, 
                                    params.gene_list ? file(params.gene_list) : empty_gene_list, 
                                    phased_bam_ch)
            } else {
                println "Phased BAM file not provided, cannot plot haplotagged DMRs"
                exit 1
            }
        } else if (bam_files_provided) {
            plot_dmrs(annotated_dmrs_ch, 
                     modbamtools_gtf, 
                     params.gene_list ? file(params.gene_list) : empty_gene_list, 
                     input_bam_ch1, 
                     input_bam_ch2)
        } else {
            println "BAM files not provided, cannot plot DMRs"
            exit 1
        }
    } else {
        if (params.'5mC') {
            println "Analysing 5mC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling_5mC(group1_beds, group2_beds)
            } else {
                prep_bedfile_results1 = prep_bedmethyl_5mC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_5mC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
            }
        } else if (params.'6mA') {
            println "Analysing 6mA methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling_6mA(group1_beds, group2_beds)
            } else {  
                prep_bedfile_results1 = prep_bedmethyl_6mA_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_6mA_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)    
            }
        } else if (params.'4mC') {
            println "Analysing 4mC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling_5mC(group1_beds, group2_beds)
            } else {  
                prep_bedfile_results1 = prep_bedmethyl_4mC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_4mC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)       
            }
        } else if (params.'5hmC') {
            println "Analysing 5hmC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling_5mC(group1_beds, group2_beds)
            } else {
                prep_bedfile_results1 = prep_bedmethyl_5hmC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_5hmC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
            }
        } else if (params.'phased_mC') {
            println "Analysing haplotagged mC methylation..."
            prep_bedfile_results1 = prep_bedmethyl_mC_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_mC_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else if (params.'phased_hmC') {
            println "Analysing haplotagged 5hmC methylation..."
            prep_bedfile_results1 = prep_bedmethyl_hmC_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_hmC_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else if (params.'phased_mA') {
            println "Analysing haplotagged 6mA methylation..."
            prep_bedfile_results1 = prep_bedmethyl_6mA_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_6mA_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else {
            println "Provide a valid methylation type flag '--5mC, --5hmC, --6mA, --4mC, --phased_mC, --phased_hmC, --phased_mA' and retry, exiting."
            exit 1
        }

        annotated_dmr_beds = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
        
        if (params.'phased_mC' || params.'phased_mA') {
            if (phased_bam_provided) {
                plot_haplotagged_dmrs(annotated_dmr_beds, modbamtools_gtf, params.gene_list ? file(params.gene_list) : empty_gene_list, phased_bam_ch)
            } else {
                println "Phased BAM file not provided, skipping plotting haplotagged DMRs"
            }
        } else if (bam_files_provided) {
            plot_dmrs(annotated_dmr_beds, modbamtools_gtf, params.gene_list ? file(params.gene_list) : empty_gene_list, input_bam_ch1, input_bam_ch2)
        } else {
            println "No BAM files provided, skipping plotting DMRs"
        }
    }
}