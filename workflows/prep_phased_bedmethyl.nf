nextflow.enable.dsl = 2

// Processing haplotagged methylcytosine (mC) file
process prep_phased_bedmethyl_mC_HP {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_mC_phased_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_mC_phased_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "m" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}

// Processing haplotagged 5-hydroxymethylcytosine (5hmC) file
process prep_phased_bedmethyl_hmC_HP {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed

  output:
    path("${bed.baseName}_hmC_phased_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_hmC_phased_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "h" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}

// Processing haplotagged 6mA methylbed file
process prep_phased_bedmethyl_mA_HP {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/modified_beds/", mode: 'copy'
  cpus 2
  memory '8 GB'
  time = 1.hours

  input:
    path bed
  output:
    path("${bed.baseName}_mA_phased_modified.bed"), emit: modified_bed
  script:
  def out_fname = "${bed.baseName}_mA_phased_modified.bed"
  """
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "a" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > ${out_fname}
  """
}
