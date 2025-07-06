nextflow.enable.dsl = 2

process dmr_calling {
  label 'ont_methyl_analysis'
  cpus 4
  memory '4 GB'
  maxForks 10  // Limit parallel chromosomes
  
  input:
    tuple val(chr), path(bed1), path(bed2)
    
  output:
    path("dmrs_${chr}.bed"), emit: chr_dmrs
    path("dmr_status_${chr}.log"), emit: status_log
    path("debug_${chr}"), emit: debug_output optional true
    
  script:
  """
  # Create debug directory for this chromosome
  mkdir -p debug_${chr}
  
  # Header for bed files
  echo -e "chr\\tpos\\tN\\tX" > header.txt
  
  # Prepare input files with header
  cat header.txt ${bed1} > ${bed1.baseName}_prepped.bed
  cat header.txt ${bed2} > ${bed2.baseName}_prepped.bed
  
  # Count sites
  n_sites_1=\$(tail -n +2 ${bed1.baseName}_prepped.bed | wc -l)
  n_sites_2=\$(tail -n +2 ${bed2.baseName}_prepped.bed | wc -l)
  
  echo "Chromosome ${chr}: Sample1 has \${n_sites_1} sites, Sample2 has \${n_sites_2} sites" > dmr_status_${chr}.log
  
  # Skip if too few sites
  if [ \${n_sites_1} -lt 10 ] || [ \${n_sites_2} -lt 10 ]; then
    echo "Insufficient sites for DMR analysis on chromosome ${chr}" >> dmr_status_${chr}.log
    echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
    exit 0
  fi
  
  # Create R script for this chromosome
  cat <<EOF > dmr_analysis_${chr}.R
  library(DSS)
  require(bsseq)

  dat1 = read.table("${bed1.baseName}_prepped.bed", header=TRUE)
  dat2 = read.table("${bed2.baseName}_prepped.bed", header=TRUE)
  
  BSobj = makeBSseqData(list(dat1, dat2), c("C1","C2"))
  dmlTest = DMLtest(BSobj, group1="C1", group2="C2", smoothing=TRUE, ncores=${task.cpus})
  
  dmrs = callDMR(dmlTest, 
                 delta=0.10, 
                 p.threshold=0.01, 
                 minlen=100, 
                 minCG=10, 
                 dis.merge=100, 
                 pct.sig=0.5)
  write.table(dmrs, file="dmrs_${chr}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
EOF

  # Run R script and capture output
  Rscript dmr_analysis_${chr}.R > r_output_${chr}.log 2>&1
  
  # Check results
  no_dmr_found=false
  if grep -q "No DMR found! Please use less stringent criteria." r_output_${chr}.log; then
    echo "No DMR found for chromosome ${chr}! Please use less stringent criteria." >> dmr_status_${chr}.log
    no_dmr_found=true
  fi
  
  # Process results
  if [ "\$no_dmr_found" = false ]; then
    if [[ -f dmrs_${chr}.tsv && -s dmrs_${chr}.tsv ]]; then
      # Add header if file has content
      echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
      tail -n +2 dmrs_${chr}.tsv >> dmrs_${chr}.bed
      n_dmrs=\$(tail -n +2 dmrs_${chr}.bed | wc -l)
      echo "Found \${n_dmrs} DMRs for chromosome ${chr}" >> dmr_status_${chr}.log
    else
      echo "No significant DMRs found for chromosome ${chr} (empty result)" >> dmr_status_${chr}.log
      echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
    fi
  else
    # Still create output file even if warning was present
    echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
    if [[ -f dmrs_${chr}.tsv && -s dmrs_${chr}.tsv ]]; then
      tail -n +2 dmrs_${chr}.tsv >> dmrs_${chr}.bed
    fi
  fi
  
  # Move debug files
  mv dmr_analysis_${chr}.R r_output_${chr}.log debug_${chr}/
  
  # Clean up
  rm -f ${bed1.baseName}_prepped.bed ${bed2.baseName}_prepped.bed header.txt dmrs_${chr}.tsv
  """
}

process split_group_beds_by_chr {
  label 'ont_methyl_analysis'
  cpus 1
  memory '2 GB'
  
  input:
    path group_beds
    val group_name
    val meth_char
    val min_cov
    
  output:
    path "chr*_${group_name}_*.bed", emit: chr_beds
    path "${group_name}_split_summary.log", emit: split_log
    
  script:
  """
  mkdir -p split_beds
  echo "Splitting ${group_name} bed files by chromosome" > ${group_name}_split_summary.log
  
  # Get all unique chromosomes across all files
  for bed_file in ${group_beds}; do
    sed 's/ /\\t/g' \${bed_file} | \
    awk -F'\\t' -v mc="${meth_char}" -v ms="${min_cov}" 'BEGIN{OFS="\\t"} \$4 == mc && \$5 >= ms {print \$1}' | \
    sort -u >> all_chromosomes.tmp
  done
  
  # Get unique chromosome list - ONLY chr1-22 and chrX
  sort -u all_chromosomes.tmp | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X)\$' > unique_chromosomes.txt
  n_chr=\$(wc -l < unique_chromosomes.txt)
  echo "Found \${n_chr} chromosomes in ${group_name} (considering only chr1-22 and chrX)" >> ${group_name}_split_summary.log
  
  # Show which chromosomes were excluded
  n_excluded=\$(sort -u all_chromosomes.tmp | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X)\$' | wc -l)
  if [ \${n_excluded} -gt 0 ]; then
    echo "Excluded \${n_excluded} non-standard chromosomes:" >> ${group_name}_split_summary.log
    sort -u all_chromosomes.tmp | grep -vE '^chr([1-9]|1[0-9]|2[0-2]|X)\$' | while read excluded_chr; do
      echo "  - \${excluded_chr}" >> ${group_name}_split_summary.log
    done
  fi
  
  # Split each file by chromosome
  while read chr; do
    for bed_file in ${group_beds}; do
      base_name=\$(basename \${bed_file} .bed)
      
      # Extract chromosome-specific data
      sed 's/ /\\t/g' \${bed_file} | \
      awk -F'\\t' -v mc="${meth_char}" -v ms="${min_cov}" -v chr="\${chr}" \
      'BEGIN{OFS="\\t"} \$1 == chr && \$4 == mc && \$5 >= ms {print \$1, \$3, \$5, \$12}' \
      > "\${chr}_${group_name}_\${base_name}.bed"
      
      # Only keep if file has content
      if [ ! -s "\${chr}_${group_name}_\${base_name}.bed" ]; then
        rm "\${chr}_${group_name}_\${base_name}.bed"
      else
        n_sites=\$(wc -l < "\${chr}_${group_name}_\${base_name}.bed")
        echo "  \${chr} - \${base_name}: \${n_sites} sites" >> ${group_name}_split_summary.log
      fi
    done
  done < unique_chromosomes.txt
  
  rm -f all_chromosomes.tmp unique_chromosomes.txt
  """
}

process group_dmr_calling {
  label 'ont_methyl_analysis'
  cpus 8
  memory '8 GB'
  maxForks 5  // Limit parallel chromosomes
  
  input:
    tuple val(chr), path(group1_beds), path(group2_beds)
    val meth_label
    
  output:
    path("dmrs_${chr}.bed"), emit: chr_dmrs
    path("dmr_status_${chr}.log"), emit: status_log
    path("debug_${chr}"), emit: debug_output optional true
    
  script:
  """
  mkdir -p g1_chr_data g2_chr_data debug_${chr} R_scripts_for_debug
  
  echo -e "chr\\tpos\\tN\\tX" > header.tsv
  echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
  echo "Processing chromosome ${chr} for ${meth_label}" > dmr_status_${chr}.log
  
  # Count samples in each group
  n_group1=\$(ls -1 ${group1_beds} 2>/dev/null | wc -l)
  n_group2=\$(ls -1 ${group2_beds} 2>/dev/null | wc -l)
  
  echo "Group 1: \${n_group1} samples, Group 2: \${n_group2} samples" >> dmr_status_${chr}.log
  
  # Prepare group1 data
  group1_has_data=false
  for bed_file in ${group1_beds}; do
    if [ -s "\${bed_file}" ]; then
      base_name=\$(basename \${bed_file} .bed | sed 's/^${chr}_group1_//')
      cat header.tsv \${bed_file} > "g1_chr_data/\${base_name}.prepped.bed"
      group1_has_data=true
    fi
  done
  
  # Prepare group2 data
  group2_has_data=false
  for bed_file in ${group2_beds}; do
    if [ -s "\${bed_file}" ]; then
      base_name=\$(basename \${bed_file} .bed | sed 's/^${chr}_group2_//')
      cat header.tsv \${bed_file} > "g2_chr_data/\${base_name}.prepped.bed"
      group2_has_data=true
    fi
  done
  
  # Only run if both groups have data
  if [ "\$group1_has_data" = true ] && [ "\$group2_has_data" = true ]; then
    cat <<'EOF' > dmr_analysis_${chr}.R
library(DSS)
require(bsseq)

  dat1_files = Sys.glob("g1_chr_data/*.prepped.bed")
  dat2_files = Sys.glob("g2_chr_data/*.prepped.bed")

  dat1_list = lapply(dat1_files, read.table, header=TRUE)
  dat2_list = lapply(dat2_files, read.table, header=TRUE)
  
  BSobj = makeBSseqData(c(dat1_list, dat2_list), 
                        c(paste0("G1_", 1:length(dat1_list)), 
                          paste0("G2_", 1:length(dat2_list))))
  
  dmlTest = DMLtest(BSobj, 
                    group1=paste0("G1_", 1:length(dat1_list)), 
                    group2=paste0("G2_", 1:length(dat2_list)), 
                    smoothing=TRUE, 
                    ncores=${task.cpus})
  
  dmrs = callDMR(dmlTest, 
                 delta=0.10, 
                 p.threshold=0.01, 
                 minlen=100, 
                 minCG=10, 
                 dis.merge=100, 
                 pct.sig=0.5)
  
  write.table(dmrs, file="dmrs_${chr}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
EOF
    
    # Run R script
    Rscript dmr_analysis_${chr}.R > r_output_${chr}.log 2>&1
    
    # Check results
    if grep -q "No DMR found! Please use less stringent criteria." r_output_${chr}.log; then
      echo "No DMR found for chromosome ${chr}! Please use less stringent criteria." >> dmr_status_${chr}.log
    fi
    
    if [[ -f dmrs_${chr}.tsv && -s dmrs_${chr}.tsv ]]; then
      tail -n +2 dmrs_${chr}.tsv >> dmrs_${chr}.bed
      n_dmrs=\$(tail -n +2 dmrs_${chr}.bed | wc -l)
      echo "Found \${n_dmrs} DMRs for chromosome ${chr}" >> dmr_status_${chr}.log
    else
      echo "No DMRs found for chromosome ${chr}" >> dmr_status_${chr}.log
    fi
    
    # Save debug info
    mv dmr_analysis_${chr}.R r_output_${chr}.log debug_${chr}/
    
  else
    echo "Skipping chromosome ${chr} - insufficient data in one or both groups" >> dmr_status_${chr}.log
  fi
  
  # Cleanup
  rm -rf g1_chr_data g2_chr_data header.tsv dmrs_${chr}.tsv
  """
}