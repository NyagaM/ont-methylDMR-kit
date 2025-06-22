nextflow.enable.dsl = 2

// Process 5mC methylbed file
process prep_bedmethyl_5mC {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours
  input:
    path bed
  output:
    path("${bed.baseName}_5mC_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_5mC_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "m" && \$5 >= 2 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}

// Process 5-hydroxymethylcytosine (5hmC) methylbed file
process prep_bedmethyl_5hmC {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours
  input:
    path bed
  output:
    path("${bed.baseName}_5hmC_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_5hmC_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "h" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}

// Process 6mA methylbed file
process prep_bedmethyl_6mA {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours
  input:
    path bed
  output:
    path("${bed.baseName}_6mA_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_6mA_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "a" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}

// Process 4mC methylbed file
process prep_bedmethyl_4mC {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours
  input:
    path bed
  output:
    path("${bed.baseName}_4mC_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_4mC_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "m" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}
