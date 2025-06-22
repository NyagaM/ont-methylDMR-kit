nextflow.enable.dsl = 2

// Evaluate significant DMRs using DSS
process dmr_calling {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/dmrs", mode: 'copy'
  cpus 24
  memory '48 GB'
  //time 2.hours

  input:
    path bed1
    path bed2

  output:
    path("dmrs_table.bed"), emit: dmr_beds
    path("dmr_status.log"), emit: status_log
    path("debug_files"), emit: debug_output_dir optional true
    
  script:
  """
  # Create debug directory
  mkdir -p debug_files R_scripts_for_debug
  
  # Header for bed files
  echo -e "chr\\tpos\\tN\\tX" > header.txt
  # Initialize an empty file for merged results
  echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed
  touch dmr_status.log
  
  # Determine available chromosomes in the input BED files
  # Use process substitution and tr to handle list for loop
  available_chromosomes=\$( (awk '{print \$1}' ${bed1} ${bed2} 2>/dev/null || true) | sort -k1,1V | uniq | grep -E '^(chr)?([0-9]+|[XYM])\$' | tr '\\n' ' ')
  
  if [ -z "\$available_chromosomes" ]; then
    echo "No valid chromosomes found in the input BED files." > dmr_status.log
    # dmrs_table.bed is already created with header
    exit 0
  fi
  echo -e "Found and processing chromosomes:\n\$(echo \$available_chromosomes | tr ' ' '\\n')" > dmr_status.log
  
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
dmlTest = DMLtest(BSobj, group1=("C1"), group2=("C2"), smoothing=TRUE, ncores=${task.cpus})
dmrs = callDMR(dmlTest, 
               delta=0.10, 
               p.threshold=0.01, 
               minlen=100, 
               minCG=10, 
               dis.merge=100, 
               pct.sig=0.5)
write.table(dmrs, file="dmrs_table_\${chr}.txt", row.names=FALSE, quote=FALSE, sep="\\t")
EOF

    # Run R script and capture output
    Rscript ont-methyl-kit_\${chr}.R > r_output_\${chr}.log 2>&1
    
    # Check if the R script found no DMRs based on the specific message
    no_dmr_found=false
    if grep -q "No DMR found! Please use less stringent criteria." r_output_\${chr}.log; then
      echo "No DMR found for chromosome \${chr}! Please use less stringent criteria." >> dmr_status.log
      no_dmr_found=true
    fi
    
    # Move R script and output log to debug directory
    mv ont-methyl-kit_\${chr}.R R_scripts_for_debug/
    mv r_output_\${chr}.log R_scripts_for_debug/

    # Only append success message if no DMR warning was found
    if [ "\$no_dmr_found" = false ]; then
      # Append the results to the merged file
      if [[ -f dmrs_table_\${chr}.txt && -s dmrs_table_\${chr}.txt ]]; then
        tail -n +2 dmrs_table_\${chr}.txt >> dmrs_table.bed
        echo "DMRs found and added for chromosome \${chr}." >> dmr_status.log
      elif [[ -f dmrs_table_\${chr}.txt ]]; then
        echo "No significant DMRs found for chromosome \${chr} (empty result file)." >> dmr_status.log
      else
        echo "Warning: dmrs_table_\${chr}.txt not found for chromosome \${chr}, skipping." >> dmr_status.log
      fi
    else
      # Still try to append any results even if warning was present
      if [[ -f dmrs_table_\${chr}.txt && -s dmrs_table_\${chr}.txt ]]; then
        tail -n +2 dmrs_table_\${chr}.txt >> dmrs_table.bed
      fi
    fi
    
    # Clean up per-chr files
    rm -f dmrs_table_\${chr}.txt ${bed1.baseName}_\${chr}.bed ${bed2.baseName}_\${chr}.bed
  done
  
  # Move debug files
  mv R_scripts_for_debug debug_files/
  """
}

// Evaluate significant group DMRs using DSS (consolidated process)
process group_dmr_calling {
  label 'ont_methyl_analysis'
  publishDir "${params.output_dir}/group_dmrs", mode: 'copy'
  cpus 36
  memory '48 GB'
  //time 6.hours // Adjust as needed

  input:
    path group1_beds
    path group2_beds
    val meth_char          // e.g., 'm' or 'a'
    val min_cov            // e.g., 5
    val meth_label         // e.g., "5mC_group" or "6mA_group" for dir naming

  output:
    path "dmrs_table.bed", emit: dmr_beds
    path "debug_files", emit: debug_output_dir // Directory containing debug files

  script:
  """
    mkdir -p debug_files R_scripts_for_debug

    echo -e "chr\\tpos\\tN\\tX" > header.tsv
    echo -e "chr\\tstart\\tend\\tlength\\tnSites\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_table.bed
    touch dmr_status.log

    # Create temporary directories for initial filtering
    mkdir -p group1_filtered_beds group2_filtered_beds

    # Initial filtering of group1 bed files
    for bed_file in ${group1_beds}; do
      base_name=\$(basename \${bed_file} .bed)
      sed 's/ /\\t/g' \${bed_file} | \
      awk -F'\\t' -v mc="${meth_char}" -v ms="${min_cov}" 'BEGIN{OFS="\\t"} \$4 == mc && \$5 >= ms {print \$1, \$3, \$5, \$12}' \
      > "group1_filtered_beds/\${base_name}.filtered.bed"
    done

    # Initial filtering of group2 bed files
    for bed_file in ${group2_beds}; do
      base_name=\$(basename \${bed_file} .bed)
      sed 's/ /\\t/g' \${bed_file} | \
      awk -F'\\t' -v mc="${meth_char}" -v ms="${min_cov}" 'BEGIN{OFS="\\t"} \$4 == mc && \$5 >= ms {print \$1, \$3, \$5, \$12}' \
      > "group2_filtered_beds/\${base_name}.filtered.bed"
    done

    # Determine available chromosomes from all filtered files
    all_chroms_file="all_chroms.tmp.txt"
    cat group1_filtered_beds/*.filtered.bed group2_filtered_beds/*.filtered.bed 2>/dev/null | awk '{print \$1}' > \$all_chroms_file

    if [ -s "\$all_chroms_file" ]; then
      available_chromosomes=\$(cat \$all_chroms_file | sort -k1,1V | uniq | grep -E '^(chr)?([0-9]+|[XYM])\$' | tr '\\n' ' ')
    else
      available_chromosomes=""
    fi
    rm -f \$all_chroms_file

    if [ -z "\$available_chromosomes" ]; then
      echo "No valid chromosomes found in input BED files after filtering for ${meth_label}." > dmr_status.log
      mv dmr_status.log debug_files/
      exit 0 # dmrs_table.bed already created with header
    fi
    echo -e "Found and processing chromosomes for ${meth_label}:\n\$(echo \$available_chromosomes | tr ' ' '\\n')" >> dmr_status.log

    # Process each chromosome from start to finish before moving to the next
    for chr_val in \$available_chromosomes; do
      echo "Processing chromosome: \${chr_val} for ${meth_label}" >&2
      
      # Create temporary directories for this chromosome
      mkdir -p "g1_chr_data_\${chr_val}" "g2_chr_data_\${chr_val}"
      
      # Prepare group1 data for the current chromosome
      group1_has_data=false
      for filtered_bed in group1_filtered_beds/*.filtered.bed; do
        base_name=\$(basename \$filtered_bed .filtered.bed)
        if [ -s "\$filtered_bed" ]; then
          awk -v c="\$chr_val" 'BEGIN{OFS="\\t"} \$1 == c {print \$0}' \$filtered_bed > "g1_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed"
          if [ -s "g1_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed" ]; then
            cat header.tsv "g1_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed" > "g1_chr_data_\${chr_val}/\${base_name}.prepped.bed"
            group1_has_data=true
          fi
          rm -f "g1_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed"
        fi
      done

      # Prepare group2 data for the current chromosome
      group2_has_data=false
      for filtered_bed in group2_filtered_beds/*.filtered.bed; do
        base_name=\$(basename \$filtered_bed .filtered.bed)
        if [ -s "\$filtered_bed" ]; then
          awk -v c="\$chr_val" 'BEGIN{OFS="\\t"} \$1 == c {print \$0}' \$filtered_bed > "g2_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed"
          if [ -s "g2_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed" ]; then
            cat header.tsv "g2_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed" > "g2_chr_data_\${chr_val}/\${base_name}.prepped.bed"
            group2_has_data=true
          fi
          rm -f "g2_chr_data_\${chr_val}/\${base_name}.chr_specific.tmp.bed"
        fi
      done
      
      # Only run R script if both groups have data for this chromosome
      if [ "\$group1_has_data" = true ] && [ "\$group2_has_data" = true ]; then
        # Create and run R script for this chromosome
        cat <<EOF > ont-methyl-kit_grouped_\${chr_val}.R
library(DSS)
require(bsseq)

dat1_files = Sys.glob("g1_chr_data_\${chr_val}/*.prepped.bed")
dat2_files = Sys.glob("g2_chr_data_\${chr_val}/*.prepped.bed")

if (length(dat1_files) == 0 || length(dat2_files) == 0) {
  warning(paste0("No prepped files found for chromosome \${chr_val} in one or both groups. Skipping this chromosome."))
} else {
  writeLines(dat1_files, con="R_scripts_for_debug/dat1_files_\${chr_val}.txt")
  writeLines(dat2_files, con="R_scripts_for_debug/dat2_files_\${chr_val}.txt")
  
  dat1_list = lapply(dat1_files, read.table, header=TRUE)
  dat2_list = lapply(dat2_files, read.table, header=TRUE)
  
  if (length(dat1_list) == 0 || length(dat2_list) == 0) {
    warning(paste0("No data read for chromosome \${chr_val} in one or both groups after attempting to read files. Skipping this chromosome."))
  } else {
    BSobj = makeBSseqData(c(dat1_list, dat2_list), 
                          c(paste0("C", 1:length(dat1_list)), paste0("N", 1:length(dat2_list))))
    
    dmlTest = DMLtest(BSobj, 
                      group1=paste0("C", 1:length(dat1_list)), 
                      group2=paste0("N", 1:length(dat2_list)), 
                      smoothing=TRUE, 
                      ncores=${task.cpus})
    
    dmrs = callDMR(dmlTest, 
                   delta=0.10, 
                   p.threshold=0.01, 
                   minlen=100, 
                   minCG=10, 
                   dis.merge=100, 
                   pct.sig=0.5)
    write.table(dmrs, file="dmrs_table_\${chr_val}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
  }
}
EOF
        
        # Run R script for this chromosome and capture output
        Rscript ont-methyl-kit_grouped_\${chr_val}.R > r_output_\${chr_val}.log 2>&1
        
        # Check if the R script found no DMRs based on the specific message
        no_dmr_found=false
        if grep -q "No DMR found! Please use less stringent criteria." r_output_\${chr_val}.log; then
          echo "No DMR found for chromosome \${chr_val}! Please use less stringent criteria." >> dmr_status.log
          no_dmr_found=true
        fi
        
        # Move R script and output log to debug directory
        mv ont-methyl-kit_grouped_\${chr_val}.R R_scripts_for_debug/
        mv r_output_\${chr_val}.log R_scripts_for_debug/
        
        # Only append success message if no DMR warning was found
        if [ "\$no_dmr_found" = false ]; then
          # Append the results to the merged file
          if [[ -f dmrs_table_\${chr_val}.tsv && -s dmrs_table_\${chr_val}.tsv ]]; then
            tail -n +2 dmrs_table_\${chr_val}.tsv >> dmrs_table.bed
            echo "DMRs found and added for chromosome \${chr_val}." >> dmr_status.log
          elif [[ -f dmrs_table_\${chr_val}.tsv ]]; then
            echo "No significant DMRs found for chromosome \${chr_val} (empty result file)." >> dmr_status.log
          else
            echo "Warning: dmrs_table_\${chr_val}.tsv not found for chromosome \${chr_val}, skipping append." >> dmr_status.log
          fi
        else
          # Still try to append any results even if warning was present
          if [[ -f dmrs_table_\${chr_val}.tsv && -s dmrs_table_\${chr_val}.tsv ]]; then
            tail -n +2 dmrs_table_\${chr_val}.tsv >> dmrs_table.bed
          fi
        fi
        
        # Clean up per-chromosome files immediately
        rm -f dmrs_table_\${chr_val}.tsv
      else
        echo "Skipping chromosome \${chr_val} - insufficient data in one or both groups." >> dmr_status.log
      fi
      
      # Clean up chromosome-specific directories immediately after processing
      rm -rf "g1_chr_data_\${chr_val}" "g2_chr_data_\${chr_val}"
    done

    # Move status log and R scripts to debug directory
    mv dmr_status.log R_scripts_for_debug/ 
    mv R_scripts_for_debug debug_files/
    
    # Final cleanup of filtered beds and header
    rm -rf group1_filtered_beds group2_filtered_beds header.tsv
  """     
}