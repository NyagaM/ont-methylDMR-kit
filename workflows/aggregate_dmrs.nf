process aggregate_dmrs {
  label 'ont_methyl_analysis'
  label 'process_low'
  publishDir "${params.output_dir}/dmrs", mode: 'copy'
  
  input:
    path chr_dmr_files
    path status_logs
    path debug_dirs
    
  output:
    path "dmrs_table.bed", emit: dmr_beds
    path "dmr_status.log", emit: status_log
    path "dmr_summary_stats.tsv", emit: summary_stats
    path "debug_files", emit: debug_output_dir
    
  script:
  """
  # Create combined debug directory
  mkdir -p debug_files
  
  # Move all chromosome debug directories
  for debug_dir in debug_chr*; do
    if [ -d "\${debug_dir}" ]; then
      mv \${debug_dir}/* debug_files/ 2>/dev/null || true
    fi
  done
  
  # Create header
  echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed
  
  # Combine all status logs
  echo "DMR Analysis Summary" > dmr_status.log
  echo "===================" >> dmr_status.log
  cat dmr_status_chr*.log >> dmr_status.log
  
  # Create summary statistics
  echo -e "chromosome\\tn_dmrs\\ttotal_length\\tmean_diff_methy\\tmax_diff_methy" > dmr_summary_stats.tsv
  
  # Process each chromosome DMR file
  for dmr_file in dmrs_chr*.bed; do
    if [[ -f "\${dmr_file}" && -s "\${dmr_file}" ]]; then
      # Extract chromosome name
      chr=\$(basename \${dmr_file} | sed 's/dmrs_//; s/.bed//')
      
      # Count DMRs (excluding header)
      n_dmrs=\$(tail -n +2 \${dmr_file} | wc -l)
      
      if [ \${n_dmrs} -gt 0 ]; then
        # Append to combined file (skip header)
        tail -n +2 \${dmr_file} >> dmrs_table.bed
        
        # Calculate statistics
        total_length=\$(tail -n +2 \${dmr_file} | awk '{sum+=\$4} END {print sum ? sum : 0}')
        mean_diff=\$(tail -n +2 \${dmr_file} | awk '{sum+=sqrt(\$8*\$8); n++} END {print n ? sum/n : 0}')
        max_diff=\$(tail -n +2 \${dmr_file} | awk 'BEGIN{max=0} {if(sqrt(\$8*\$8)>max) max=sqrt(\$8*\$8)} END {print max}')
        
        echo -e "\${chr}\\t\${n_dmrs}\\t\${total_length}\\t\${mean_diff}\\t\${max_diff}" >> dmr_summary_stats.tsv
      else
        echo -e "\${chr}\\t0\\t0\\t0\\t0" >> dmr_summary_stats.tsv
      fi
    fi
  done
  
  # Sort DMRs by chromosome and position
  if [ \$(tail -n +2 dmrs_table.bed | wc -l) -gt 0 ]; then
    head -1 dmrs_table.bed > dmrs_table_sorted.bed
    tail -n +2 dmrs_table.bed | sort -k1,1V -k2,2n >> dmrs_table_sorted.bed
    mv dmrs_table_sorted.bed dmrs_table.bed
  fi
  
  # Final summary
  total_dmrs=\$(tail -n +2 dmrs_table.bed | wc -l)
  echo "" >> dmr_status.log
  echo "===================" >> dmr_status.log
  echo "Total DMRs across all chromosomes: \${total_dmrs}" >> dmr_status.log
  
  # Create a visual summary
  echo "" >> dmr_status.log
  echo "DMRs per chromosome:" >> dmr_status.log
  tail -n +2 dmr_summary_stats.tsv | awk '{printf "  %s: %d DMRs\\n", \$1, \$2}' >> dmr_status.log
  """
}