process prep_phased_bedmethyl_mC_HP {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/modified_beds/mC_phased", mode: 'copy'

  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_mC_phased_modified.bed", emit: chr_beds
    path "${bed.baseName}_mC_phased_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing phased mC bedmethyl file: ${bed}" > ${bed.baseName}_mC_phased_prep_summary.log
  echo "Filtering for mC sites with score >= 5" >> ${bed.baseName}_mC_phased_prep_summary.log
  
  # First, create a filtered bed with only mC sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "m" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_mC_phased_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_mC_phased_temp.bed)
  echo "Total phased mC sites found: \${total_sites}" >> ${bed.baseName}_mC_phased_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_mC_phased_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_mC_phased_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_mC_phased_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_mC_phased_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_mC_phased_temp.bed > "\${chr}_${bed.baseName}_mC_phased_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_mC_phased_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_mC_phased_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_mC_phased_temp.bed chromosomes.txt
  """
}

process prep_phased_bedmethyl_hmC_HP {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/modified_beds/hmC_phased", mode: 'copy'
  
  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_hmC_phased_modified.bed", emit: chr_beds
    path "${bed.baseName}_hmC_phased_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing phased hmC bedmethyl file: ${bed}" > ${bed.baseName}_hmC_phased_prep_summary.log
  echo "Filtering for hmC sites with score >= 5" >> ${bed.baseName}_hmC_phased_prep_summary.log
  
  # First, create a filtered bed with only hmC sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "h" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_hmC_phased_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_hmC_phased_temp.bed)
  echo "Total phased hmC sites found: \${total_sites}" >> ${bed.baseName}_hmC_phased_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_hmC_phased_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_hmC_phased_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_hmC_phased_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_hmC_phased_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_hmC_phased_temp.bed > "\${chr}_${bed.baseName}_hmC_phased_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_hmC_phased_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_hmC_phased_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_hmC_phased_temp.bed chromosomes.txt
  """
}

process prep_phased_bedmethyl_mA_HP {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/modified_beds/mA_phased", mode: 'copy'
  
  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_mA_phased_modified.bed", emit: chr_beds
    path "${bed.baseName}_mA_phased_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing phased mA bedmethyl file: ${bed}" > ${bed.baseName}_mA_phased_prep_summary.log
  echo "Filtering for mA sites with score >= 5" >> ${bed.baseName}_mA_phased_prep_summary.log
  
  # First, create a filtered bed with only mA sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "a" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_mA_phased_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_mA_phased_temp.bed)
  echo "Total phased mA sites found: \${total_sites}" >> ${bed.baseName}_mA_phased_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_mA_phased_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_mA_phased_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_mA_phased_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_mA_phased_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_mA_phased_temp.bed > "\${chr}_${bed.baseName}_mA_phased_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_mA_phased_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_mA_phased_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_mA_phased_temp.bed chromosomes.txt
  """
}