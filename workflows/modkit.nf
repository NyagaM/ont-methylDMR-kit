process modkit_mC {
    label "modkit"
    label "process_high"
    publishDir "${params.output_dir}/modkit/5mC", mode: 'copy'
   
    input:
        tuple val(name), path(bam), path(bai)
        tuple val(ref_name), path(reference), path(fai)
        
    output:
        tuple val(name), path("${name}*.bed"), emit: modkit
        tuple val(name), path("*_${name}.log"), emit: modkit_log

    script:
    """
    # compute sample probabilities
    probs=\$( modkit sample-probs ${bam} \\
        -p 0.1 \\
        --interval-size 5000000 \\
        --only-mapped \\
        --threads ${task.cpus} \\
        2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )

    # run modkit pileup for 5mC
    modkit pileup \\
        --ref ${reference} \\
        --interval-size 1000000 \\
        --log-filepath modkit_${name}.log \\
        ${probs} \\
        --combine-strands \\
        --cpg \\
        --threads ${task.cpus} \\
        ${bam} \\
        ${name}_5mC.bed
    """
}

process modkit_phased_mC {
    label "modkit"
    publishDir "${params.output_dir}/modkit/phased_5mC", mode: 'copy'
    
    input:
        tuple val(name), path(bam), path(bai)
        tuple val(ref_name), path(reference), path(fai)

    output:
        tuple val(name), path("${name}*ungrouped.bed"), emit: modkit_H0, optional: true
        tuple val(name), path("${name}*1.bed"), emit: modkit_H1, optional: true
        tuple val(name), path("${name}*2.bed"), emit: modkit_H2, optional: true
        tuple val(name), path("*_${name}.log"), emit: modkit_log

    script:
    """
    # compute sample probabilities
    probs=\$( modkit sample-probs ${bam} \\
        -p 0.1 \\
        --interval-size 5000000 \\
        --only-mapped \\
        --threads ${task.cpus} \\
        2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )

    # run modkit pileup for phased 5mC
    modkit pileup \\
        --ref ${reference} \\
        --interval-size 1000000 \\
        --log-filepath modkit_phased_${name}.log \\
        ${probs} \\
        --prefix ${name}_5mC \\
        --partition-tag HP \\
        --combine-strands \\
        --cpg \\
        --threads ${task.cpus} \\
        ${bam} \\
        ${name}

    mv ${name}/*.bed ./ # Move the output files to the current directory 
    """
}

process modkit_mA {
    label "modkit"
    label "process_high"
    publishDir "${params.output_dir}/modkit/6mA", mode: 'copy'
    
    input:
        tuple val(name), path(bam), path(bai)
        tuple val(ref_name), path(reference), path(fai)
        
    output:
        tuple val(name), path("${name}*.bed"), emit: modkit
        tuple val(name), path("*_${name}.log"), emit: modkit_log

    script:
    """
    # compute sample probabilities
    probs=\$( modkit sample-probs ${bam} \\
        -p 0.1 \\
        --interval-size 5000000 \\
        --only-mapped \\
        --threads ${task.cpus} \\
        2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )

    # run modkit pileup for 6mA
    modkit pileup \\
        --ref ${reference} \\
        --interval-size 1000000 \\
        --log-filepath modkit_6mA_${name}.log \\
        ${probs} \\
        --motif AGG 0 \\
        --threads ${task.cpus} \\
        ${bam} \\
        ${name}_6mA.bed
    """
}

process modkit_phased_mA {
    label "modkit"
    publishDir "${params.output_dir}/modkit/phased_6mA", mode: 'copy'
    
    input:
        tuple val(name), path(bam), path(bai)
        tuple val(ref_name), path(reference), path(fai)
        
    output:
        tuple val(name), path("${name}*ungrouped.bed"), emit: modkit_H0, optional: true
        tuple val(name), path("${name}*1.bed"), emit: modkit_H1, optional: true
        tuple val(name), path("${name}*2.bed"), emit: modkit_H2, optional: true
        tuple val(name), path("*_${name}.log"), emit: modkit_log

    script:
    """
    # compute sample probabilities
    probs=\$( modkit sample-probs ${bam} \\
        -p 0.1 \\
        --interval-size 5000000 \\
        --only-mapped \\
        --threads ${task.cpus} \\
        2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )

    # run modkit pileup for phased 6mA
    modkit pileup \\
        --ref ${reference} \\
        --interval-size 1000000 \\
        --log-filepath modkit_phased_6mA_${name}.log \\
        ${probs} \\
        --prefix ${name}_6mA \\
        --partition-tag HP \\
        --motif AGG 0 \\
        --threads ${task.cpus} \\
        ${bam} \\
        ${name}
    
    mv ${name}/*.bed ./ # Move the output files to the current directory 
    """
}