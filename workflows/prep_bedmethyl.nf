process prep_bedmethyl_5mC {
  label 'ont_methyl_analysis'
  label 'process_default'
  publishDir "${params.output_dir}/modified_beds/5mC", mode: 'copy'
  
  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_5mC_modified.bed", emit: chr_beds
    path "${bed.baseName}_5mC_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing 5mC bedmethyl file: ${bed}" > ${bed.baseName}_5mC_prep_summary.log
  echo "Filtering for 5mC sites with read depth >= 5" >> ${bed.baseName}_5mC_prep_summary.log
  
  # First, create a filtered bed with only 5mC sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "m" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_5mC_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_5mC_temp.bed)
  echo "Total 5mC sites found: \${total_sites}" >> ${bed.baseName}_5mC_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_5mC_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_5mC_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_5mC_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_5mC_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_5mC_temp.bed > "\${chr}_${bed.baseName}_5mC_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_5mC_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_5mC_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_5mC_temp.bed chromosomes.txt
  """
}

process prep_bedmethyl_5hmC {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/modified_beds/5hmC", mode: 'copy'
  
  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_5hmC_modified.bed", emit: chr_beds
    path "${bed.baseName}_5hmC_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing 5hmC bedmethyl file: ${bed}" > ${bed.baseName}_5hmC_prep_summary.log
  echo "Filtering for 5hmC sites with read depth >= 5" >> ${bed.baseName}_5hmC_prep_summary.log
  
  # First, create a filtered bed with only 5hmC sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "h" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_5hmC_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_5hmC_temp.bed)
  echo "Total 5hmC sites found: \${total_sites}" >> ${bed.baseName}_5hmC_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_5hmC_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_5hmC_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_5hmC_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_5hmC_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_5hmC_temp.bed > "\${chr}_${bed.baseName}_5hmC_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_5hmC_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_5hmC_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_5hmC_temp.bed chromosomes.txt
  """
}

process prep_bedmethyl_6mA {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/modified_beds/6mA", mode: 'copy'
  
  input:
    path bed
    
  output:
    path "chr*_${bed.baseName}_6mA_modified.bed", emit: chr_beds
    path "${bed.baseName}_6mA_prep_summary.log", emit: prep_log
    
  script:
  """
  echo "Preparing 6mA bedmethyl file: ${bed}" > ${bed.baseName}_6mA_prep_summary.log
  echo "Filtering for 6mA sites with read depth >= 5" >> ${bed.baseName}_6mA_prep_summary.log
  
  # First, create a filtered bed with only 6mA sites
  sed 's/ /\\t/g' ${bed} | awk -F'\\t' 'BEGIN{OFS="\\t"} \$4 == "a" && \$5 >= 5 {print \$1, \$3, \$5, \$12}' > filtered_6mA_temp.bed
  
  # Get total sites
  total_sites=\$(wc -l < filtered_6mA_temp.bed)
  echo "Total 6mA sites found: \${total_sites}" >> ${bed.baseName}_6mA_prep_summary.log
  
  # Get list of chromosomes - chr1-22, chrX and chrM
  cut -f1 filtered_6mA_temp.bed | sort -V | uniq | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' > chromosomes.txt
  n_chr=\$(wc -l < chromosomes.txt)
  echo "Found \${n_chr} chromosomes (Checked for chr1-22, chrX and chrM only)" >> ${bed.baseName}_6mA_prep_summary.log
  
  # Show excluded chromosomes
  n_excluded=\$(cut -f1 filtered_6mA_temp.bed | sort -V | uniq | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X|M)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes" >> ${bed.baseName}_6mA_prep_summary.log
  fi
  
  # Split by chromosome
  while read chr; do
    awk -v chr="\${chr}" 'BEGIN{OFS="\\t"} \$1 == chr' filtered_6mA_temp.bed > "\${chr}_${bed.baseName}_6mA_modified.bed"
    n_sites=\$(wc -l < "\${chr}_${bed.baseName}_6mA_modified.bed")
    echo "  \${chr}: \${n_sites} sites" >> ${bed.baseName}_6mA_prep_summary.log
  done < chromosomes.txt
  
  # Clean up
  rm filtered_6mA_temp.bed chromosomes.txt
  """
}
