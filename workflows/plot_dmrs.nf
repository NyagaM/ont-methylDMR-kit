nextflow.enable.dsl = 2

process plot_dmr_modbamtools {
  label 'modbamtools'
  publishDir "${params.output_dir}/dmr_plots/modbamtools", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    tuple file(gtf_file), file(gtf_tbi_file)
    path gene_list_path_input
    tuple val(name1), file(bam1), file(bai1)
    tuple val(name2), file(bam2), file(bai2)

  output:
    path("*.html"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs
    path("plot_status.log"), emit: status_log, optional: true
    
  script:
  """
  genes_of_interest_source_file="${gene_list_path_input}"
  genes_of_interest_local_copy="genes_of_interest.txt"

  if [[ -f "\${genes_of_interest_source_file}" && -s "\${genes_of_interest_source_file}" ]]; then
    cp "\${genes_of_interest_source_file}" "\${genes_of_interest_local_copy}"
  else
    touch "\${genes_of_interest_local_copy}"
  fi

  # Create an empty error log file and status log
  touch error.log
  touch plot_status.log

  # Use a temporary file to track plot count
  echo "0" > plot_count.tmp

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nSites meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s "\${genes_of_interest_local_copy}" ]] || grep -q -w "\${gene}" "\${genes_of_interest_local_copy}"; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      
      if modbamtools plot -r \${region} -g ${gtf_file} -s ${bam1.baseName},${bam2.baseName} -p \${output_prefix} ${bam1} ${bam2} -o ./; then
        # Increment plot count on success
        current_count=\$(cat plot_count.tmp)
        echo \$((current_count + 1)) > plot_count.tmp
        echo "Successfully plotted region \${region} for gene \${gene} using modbamtools" >> plot_status.log
      else
        echo "Error processing region \${region} for gene \${gene}" >> error.log
      fi
    fi
  done

  # Read the final plot count
  final_plot_count=\$(cat plot_count.tmp)
  
  echo "Total plots generated: \${final_plot_count}" >> plot_status.log
  
  if [[ \${final_plot_count} -eq 0 ]]; then
    echo "No plots generated, exiting successfully." >> error.log
    echo "No plots generated - either no DMRs matched the criteria or all plotting attempts failed." >> plot_status.log
  fi
  
  # Clean up temporary file
  rm -f plot_count.tmp
  """
}

process plot_dmr_methylartist {
  label 'methylartist'
  publishDir "${params.output_dir}/dmr_plots/methylartist", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    tuple file(gtf_file), file(gtf_tbi_file)
    path gene_list_path_input
    path reference
    tuple val(name1), file(bam1), file(bai1)
    tuple val(name2), file(bam2), file(bai2)

  output:
    path("*.png"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs, optional: true
    path("plot_status.log"), emit: status_log, optional: true
    
  script:
  """
  genes_of_interest_source_file="${gene_list_path_input}"
  genes_of_interest_local_copy="genes_of_interest.txt"

  if [[ -f "\${genes_of_interest_source_file}" && -s "\${genes_of_interest_source_file}" ]]; then
    cp "\${genes_of_interest_source_file}" "\${genes_of_interest_local_copy}"
  else
    touch "\${genes_of_interest_local_copy}"
  fi

  # Create empty log files
  touch error.log
  touch plot_status.log

  # Use a temporary file to track plot count
  echo "0" > plot_count.tmp

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nSites meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s "\${genes_of_interest_local_copy}" ]] || grep -q -w "\${gene}" "\${genes_of_interest_local_copy}"; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      
      if methylartist locus --interval \${region} --gtf ${gtf_file} --bams ${bam1},${bam2} --ref ${reference} --motif CG --mods m --outfile \${output_prefix} --labelgenes --nomask; then
        # Increment plot count on success
        current_count=\$(cat plot_count.tmp)
        echo \$((current_count + 1)) > plot_count.tmp
        echo "Successfully plotted region \${region} for gene \${gene} using methylartist" >> plot_status.log
      else
        echo "Error processing region \${region} for gene \${gene}" >> error.log
      fi
    fi
  done

  # Read the final plot count
  final_plot_count=\$(cat plot_count.tmp)
  
  echo "Total methylartist plots generated: \${final_plot_count}" >> plot_status.log
  
  if [[ \${final_plot_count} -eq 0 ]]; then
    echo "No methylartist plots generated, exiting successfully." >> error.log
    echo "No methylartist plots generated - either no DMRs matched the criteria or all plotting attempts failed." >> plot_status.log
  fi
  
  # Clean up temporary file
  rm -f plot_count.tmp
  """
}

process plot_phased_dmr_modbamtools {
  label 'modbamtools'
  publishDir "${params.output_dir}/haplotagged_dmr_plots/modbamtools", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    tuple file(gtf_file), file(gtf_tbi_file)
    path gene_list_path_input
    tuple val(name), file(bam), file(bai)

  output:
    path("*.html"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs
    path("plot_status.log"), emit: status_log, optional: true

  script:
  """
  genes_of_interest_source_file="${gene_list_path_input}"
  genes_of_interest_local_copy="genes_of_interest.txt"

  if [[ -f "\${genes_of_interest_source_file}" && -s "\${genes_of_interest_source_file}" ]]; then
    cp "\${genes_of_interest_source_file}" "\${genes_of_interest_local_copy}"
  else
    touch "\${genes_of_interest_local_copy}"
  fi

  # Create empty log files
  touch error.log
  touch plot_status.log

  # Use a temporary file to track plot count
  echo "0" > plot_count.tmp

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s "\${genes_of_interest_local_copy}" ]] || grep -q -w "\${gene}" "\${genes_of_interest_local_copy}"; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      
      if modbamtools plot -r \${region} -g ${gtf_file} -s ${bam.baseName} -p \${output_prefix} -hp ${bam} -o ./; then
        # Increment plot count on success
        current_count=\$(cat plot_count.tmp)
        echo \$((current_count + 1)) > plot_count.tmp
        echo "Successfully plotted phased region \${region} for gene \${gene} using modbamtools" >> plot_status.log
      else
        echo "Error processing phased region \${region} for gene \${gene}" >> error.log
      fi
    fi
  done

  # Read the final plot count
  final_plot_count=\$(cat plot_count.tmp)
  
  echo "Total phased plots generated: \${final_plot_count}" >> plot_status.log
  
  if [[ \${final_plot_count} -eq 0 ]]; then
    echo "No phased plots generated, exiting successfully." >> error.log
    echo "No phased plots generated - either no DMRs matched the criteria or all plotting attempts failed." >> plot_status.log
  fi
  
  # Clean up temporary file
  rm -f plot_count.tmp
  """
}

process plot_phased_dmr_methylartist {
  label 'methylartist'
  publishDir "${params.output_dir}/haplotagged_dmr_plots/methylartist", mode: 'copy'
  cpus 2
  memory '4 GB'
  time 9.hours

  input:
    path annotated_dmr_beds
    tuple file(gtf_file), file(gtf_tbi_file)
    path gene_list_path_input
    path reference
    tuple val(name), file(bam), file(bai)
    
  output:
    path("*.png"), emit: dmr_plots, optional: true
    path("error.log"), emit: dmr_logs, optional: true
    path("plot_status.log"), emit: status_log, optional: true
    
  script:
  """
  genes_of_interest_source_file="${gene_list_path_input}"
  genes_of_interest_local_copy="genes_of_interest.txt"

  if [[ -f "\${genes_of_interest_source_file}" && -s "\${genes_of_interest_source_file}" ]]; then
    cp "\${genes_of_interest_source_file}" "\${genes_of_interest_local_copy}"
  else
    touch "\${genes_of_interest_local_copy}"
  fi

  # Create empty log files
  touch error.log
  touch plot_status.log

  # Use a temporary file to track plot count
  echo "0" > plot_count.tmp

  awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$4 FS \$NF; if (!seen[key]++) print }' ${annotated_dmr_beds} | while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
    if [[ ! -s "\${genes_of_interest_local_copy}" ]] || grep -q -w "\${gene}" "\${genes_of_interest_local_copy}"; then
      region="\${chr}:\${start}-\${end}"
      region2="\${chr}_\${start}-\${end}"
      output_prefix="\${region2}_\${gene}"
      
      if methylartist locus --interval \${region} --gtf ${gtf_file} --bams ${bam} --ref ${reference} --motif CG --mods m --outfile \${output_prefix} --labelgenes --nomask --phased --ignore_ps; then
        # Increment plot count on success
        current_count=\$(cat plot_count.tmp)
        echo \$((current_count + 1)) > plot_count.tmp
        echo "Successfully plotted phased region \${region} for gene \${gene} using methylartist" >> plot_status.log
      else
        echo "Error processing region \${region} for gene \${gene}" >> error.log
      fi
    fi
  done

  # Read the final plot count
  final_plot_count=\$(cat plot_count.tmp)
  
  echo "Total methylartist phased plots generated: \${final_plot_count}" >> plot_status.log
  
  if [[ \${final_plot_count} -eq 0 ]]; then
    echo "No methylartist phased plots generated, exiting successfully." >> error.log
    echo "No methylartist phased plots generated - either no DMRs matched the criteria or all plotting attempts failed." >> plot_status.log
  fi
  
  # Clean up temporary file
  rm -f plot_count.tmp
  """
}