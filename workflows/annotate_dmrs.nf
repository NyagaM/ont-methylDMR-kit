process annotate_dmrs {
  label 'ont_methyl_analysis'
  label 'process_default'
  publishDir "${params.output_dir}/annotated_dmrs", mode: 'copy'

  input:
    path bed
    path annotationFile

  output:
    path("dmrs_table_annotated.bed"), emit: annotated_dmr_beds
    path("dmrs_table_annotated_imprinted.bed"), emit: annotated_dmr_beds_imprinted, optional: true
    path("annotation_summary.log"), emit: annotation_log

  script:
  """
  # Create header for annotated files
  echo -e 'chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat\\tannotation_chr\\tannotation_start\\tannotation_end\\tstrand\\tannotation\\tbiotype\\tgene' > annotation_header.txt
  
  # Create headerless BED file for bedtools (skip first line if it contains 'chr\tstart\tend')
  if head -n 1 ${bed} | grep -q "^chr[[:space:]]\\+start[[:space:]]\\+end"; then
    echo "Detected header in BED file, removing it for bedtools compatibility"
    tail -n +2 ${bed} > dmrs_clean.bed
  else
    echo "No header detected in BED file, using as-is"
    cp ${bed} dmrs_clean.bed
  fi
  
  # Copy annotation file to resolve any symlink issues
  cp -L ${annotationFile} annotation_clean.bed
  
  # Perform bedtools intersection with headerless BED file
  bedtools intersect -a dmrs_clean.bed -b annotation_clean.bed -wa -wb > dmrs_table_annotated.tmp
  
  # Create the full annotated DMR table with unique entries only
  (head -n 1 annotation_header.txt; sort dmrs_table_annotated.tmp | uniq) > dmrs_table_annotated.bed
  
  # Create summary log
  touch annotation_summary.log
  total_dmrs=\$(wc -l < dmrs_clean.bed)
  annotated_dmrs=\$(sort dmrs_table_annotated.tmp | uniq | wc -l)
  echo "Total DMRs: \${total_dmrs}" >> annotation_summary.log
  echo "DMR Annotations: \${annotated_dmrs}" >> annotation_summary.log
  
  # If --imprinted flag is used, create a filtered table
  if [ "${params.imprinted}" = "true" ]; then
    # Get the imprinted genes file path
    imprinted_genes_file="${workflow.projectDir}/annotations/imprinted_genes.tsv"
    
    if [ -f "\${imprinted_genes_file}" ]; then
      # Extract gene names from imprinted genes file (assuming first column after header)
      tail -n +2 "\${imprinted_genes_file}" | cut -f1 | sort | uniq > imprinted_genes_list.txt
      
      # Filter the annotated DMRs to include only those overlapping imprinted genes
      # The gene name is in the last column (field 16)
      awk -F'\\t' 'NR==FNR{genes[\$1]; next} FNR==1 || \$16 in genes' imprinted_genes_list.txt dmrs_table_annotated.bed > dmrs_table_annotated_imprinted.bed
      
      # Count imprinted DMRs (unique based on chr, start, end)
      imprinted_dmrs=\$(tail -n +2 dmrs_table_annotated_imprinted.bed | cut -f1-3 | sort | uniq | wc -l)
      echo "DMRs overlapping imprinted genes: \${imprinted_dmrs}" >> annotation_summary.log
      
      # List unique imprinted genes with DMRs
      echo -e "\\nImprinted genes with DMRs:" >> annotation_summary.log
      tail -n +2 dmrs_table_annotated_imprinted.bed | cut -f1-3,16 | sort | uniq | cut -f4 | sort | uniq | while read gene; do
        count=\$(tail -n +2 dmrs_table_annotated_imprinted.bed | cut -f1-3,16 | sort | uniq | awk -F'\\t' -v gene="\${gene}" '\$4 == gene' | wc -l)
        echo "  \${gene}: \${count} DMRs" >> annotation_summary.log
      done
      
      # Cleanup
      rm imprinted_genes_list.txt
    else
      echo "WARNING: Imprinted genes file not found at \${imprinted_genes_file}" >> annotation_summary.log
      echo "Skipping imprinted gene filtering" >> annotation_summary.log
    fi
  fi
  
  # Cleanup temporary files
  rm dmrs_table_annotated.tmp annotation_header.txt dmrs_clean.bed annotation_clean.bed
  """
}
