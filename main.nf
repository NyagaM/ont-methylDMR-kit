#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Help message
def helpMessage = """
Usage: nextflow run main.nf [options]

Options:
  --input_file1        First input bedmethyl file from sample_1 or haplotype_1 bedmethyl (required for analysis)
  --input_file2        Second input bedmethyl file from sample_2 or haplotype_2 bedmethyl (required for analysis)
  --output_dir         Output directory (required)
  --input_modBam1      First input modified BAM used to generate bedmethyl for sample_1 (optional, for plotting)
  --input_modBam2      Second input modified BAM used to generate bedmethyl for sample_2 (optional, for plotting)
  --input_group1       Directory containing bedmethyl files (with *.bed extension) for group 1 (optional, for group analysis)
  --input_group2       Directory containing bedmethyl files (with *.bed extension) for group 2 (optional, for group analysis)
  --gene_list          A list of genes to generate plots on (optional). If not provided and --imprinted is used, imprinted genes will be plotted.
  --reference          Reference genome FASTA file (required for methylartist plots)
  --imprinted          Use this flag to plot DMRs overlapping imprinted genes if --gene_list is not provided (optional, for haplotagged plots)
  --5mC                Use this flag to trigger 5mC DMR calling (optional)
  --5hmC               Use this flag to trigger 5hmC DMR calling (optional)
  --6mA                Use this flag to trigger 6mA DMR calling (optional)
  --4mC                Use this flag to trigger 4mC DMR calling (optional)
  --phased_mC          Use this flag to trigger haplotagged 5mC/4mC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_mA          Use this flag to trigger haplotagged 6mA DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_hmC         Use this flag to trigger haplotagged 5hmC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_modBam      Haplotagged modified BAM (required for plotting DMRs if --phased_mC, --phased_mA or --phased_hmC is used)
  --plots_only         Only run the plotting processes, requires --annotated_dmrs and relevant BAM files (and --reference for methylartist)
  --annotated_dmrs     Path to annotated DMR bed file (required if --plots_only is used)
  --methylartist_only  Use this flag to only generate plots using methylartist (requires --reference). Otherwise, both modbamtools and methylartist plots are generated if reference is provided.
  --help               Print this help message
"""

params.help = false
params.input_modBam1 = false
params.input_modBam2 = false
params.gene_list = false
params.reference = false
params.phased_modBam = ""
params.input_file1 = ""
params.input_file2 = ""
params.input_group1 = ""
params.input_group2 = ""
params.output_dir = ""
params.'5mC' = false
params.'5hmC' = false
params.'6mA' = false
params.'4mC' = false
params.'phased_mA' = false
params.'phased_mC' = false
params.'phased_hmC' = false
params.plots_only = false
params.annotated_dmrs = ""
params.imprinted = false
params.methylartist_only = false

// Help information
if (params.help) {
    println helpMessage
    exit 0
}

// Early validation for --reference if --methylartist_only is used
if (params.methylartist_only && !params.reference) {
    println "ERROR: --reference genome is required when --methylartist_only is specified."
    println helpMessage
    exit 1
}

// Separate validation for plots_only mode
if (params.plots_only) {
    if (!params.annotated_dmrs || !params.output_dir) {
        println "ERROR: In plots-only mode, --annotated_dmrs and --output_dir are required."
        println helpMessage
        exit 0
    }
    boolean phased_bams_for_plots = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC') && params.phased_modBam
    boolean non_phased_bams_for_plots = params.input_modBam1 && params.input_modBam2
    if (!phased_bams_for_plots && !non_phased_bams_for_plots) {
        println "ERROR: In plots-only mode, relevant BAM files (--phased_modBam or --input_modBam1/2) are required for plotting."
        println helpMessage
        exit 0
    }
    if (params.methylartist_only && !params.reference) { // This case is already caught above, but good to be explicit for plots_only
        // This specific check might be redundant due to the global check above, but doesn't hurt.
    }
} else {
    // Validation for regular analysis mode
    if ((!params.input_file1 && !params.input_group1) || 
        (!params.input_file2 && !params.input_group2) || 
        !params.output_dir) {
        println "In analysis mode, input files (--input_file1/2 or --input_group1/2) and --output_dir are required"
        println helpMessage
        exit 0
    }
    if (params.methylartist_only && !params.reference) { // Already caught by global check if --methylartist_only is true
        // This specific check might be redundant due to the global check above.
        // println "ERROR: In analysis mode, --reference is required when --methylartist_only is specified for plotting."
        // println helpMessage
        // exit 1
    }
}

// Check if output directory exists
def outputDir = file(params.output_dir)
if (outputDir.exists()) {
    println "INFO: Output directory exists: ${params.output_dir}"
} else {
    println "INFO: Output directory does not exist. Creating one: ${params.output_dir}"
    outputDir.mkdirs()
}

// Input channels
input_ch1 = params.input_file1 ? Channel.fromPath(params.input_file1, checkIfExists: true) : Channel.empty()
input_ch2 = params.input_file2 ? Channel.fromPath(params.input_file2, checkIfExists: true) : Channel.empty()
annotationFile = file("${workflow.projectDir}/annotations/gencode.v46.annotation.exon-promoters-introns.sorted.bed") // For annotate_dmrs

// Define GTF and its index for plotting processes
gencode_annotation_gtf_file = file("${workflow.projectDir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz")
gencode_annotation_gtf_tbi = file("${workflow.projectDir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz.tbi")

if (!gencode_annotation_gtf_file.exists()) {
    println "ERROR: GTF file not found: ${gencode_annotation_gtf_file}"
    exit 1
}
if (!gencode_annotation_gtf_tbi.exists()) {
    println "ERROR: GTF index file (.tbi) not found: ${gencode_annotation_gtf_tbi}. Please ensure it exists and is named correctly (e.g., your_file.gtf.gz.tbi)."
    exit 1
}
gencode_annotation_for_plotting = [gencode_annotation_gtf_file, gencode_annotation_gtf_tbi]
annotated_dmrs_ch = params.plots_only ? Channel.fromPath(params.annotated_dmrs, checkIfExists: true) : Channel.empty()

// Check for BAM files
bam_files_provided = params.input_modBam1 && params.input_modBam2
phased_bam_provided = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC') && params.phased_modBam

input_bam_ch1 = bam_files_provided ? Channel.fromPath(params.input_modBam1, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()
input_bam_ch2 = bam_files_provided ? Channel.fromPath(params.input_modBam2, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()
phased_bam_ch = phased_bam_provided ? Channel.fromPath(params.phased_modBam, checkIfExists: true).map { [it.baseName, file(it), file(it + '.bai')] } : Channel.empty()

// Create an empty gene list file if none is provided
empty_gene_list = file("${params.output_dir}/empty_gene_list.txt")
if (!params.gene_list && !params.imprinted) { // Only create if neither gene_list nor imprinted (which implies a list) is set
    empty_gene_list.text = ""
}

// Include processes for standard bedmethyl preparation
include { prep_bedmethyl_5mC  as prep_bedmethyl_5mC_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5mC  as prep_bedmethyl_5mC_2  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5hmC as prep_bedmethyl_5hmC_1 } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5hmC as prep_bedmethyl_5hmC_2 } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_6mA  as prep_bedmethyl_6mA_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_6mA  as prep_bedmethyl_6mA_2  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_4mC  as prep_bedmethyl_4mC_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_4mC  as prep_bedmethyl_4mC_2  } from './workflows/prep_bedmethyl.nf'
// Include processes for phased bedmethyl preparation
include { prep_phased_bedmethyl_mC_HP  as prep_bedmethyl_mC_HP_1  } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mC_HP  as prep_bedmethyl_mC_HP_2  } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_hmC_HP as prep_bedmethyl_hmC_HP_1 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_hmC_HP as prep_bedmethyl_hmC_HP_2 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mA_HP  as prep_bedmethyl_6mA_HP_1 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mA_HP  as prep_bedmethyl_6mA_HP_2 } from './workflows/prep_phased_bedmethyl.nf'
// Include DMR calling processes (group_dmr_calling is the new consolidated one)
include { dmr_calling; group_dmr_calling } from './workflows/dmr_calling.nf'
// Include annotating and plotting processes
include { annotate_dmrs } from './workflows/annotate_dmrs.nf'
include { plot_dmr_modbamtools; plot_dmr_methylartist; plot_phased_dmr_modbamtools; plot_phased_dmr_methylartist } from './workflows/plot_dmrs.nf'

//workflow
workflow {

    // Determine the gene list to use for plotting
    def gene_list_for_plotting
    if (params.imprinted && !params.gene_list) {
        def imprinted_file = file("${workflow.projectDir}/annotations/imprinted_genes.tsv")
        if (imprinted_file.exists()) {
            def processed_file = file("${params.output_dir}/imprinted_genes.txt")
            try {
                def lines = imprinted_file.readLines()
                if (lines.size() > 1) {
                    def genes = lines.drop(1).collect { it.split('\t')[0].trim() }.findAll { it && !it.isEmpty() }
                    processed_file.text = genes.join('\n')
                } else {
                    processed_file.text = ""
                    println "WARNING: Imprinted genes file ${imprinted_file.name} is empty or has no valid gene entries. No genes will be plotted."
                }
                gene_list_for_plotting = processed_file
                println "INFO: Using processed imprinted gene list for plotting or annotation: ${gene_list_for_plotting.name}"
            } catch (Exception e) {
                println "ERROR: Failed to process imprinted genes file ${imprinted_genes.name}: ${e.getMessage()}"
                exit 1
            }
        } else {
            println "ERROR: --imprinted flag was used, but the imprinted genes file was not found."
            exit 1
        }
    } else if (params.gene_list) {
        gene_list_for_plotting = file(params.gene_list)
        println "INFO: Using user-provided gene list for plotting: ${params.gene_list}"
    } else {
        // Ensure empty_gene_list file is created if it's going to be used
        if (!empty_gene_list.exists()) {
            empty_gene_list.text = ""
        }
        gene_list_for_plotting = empty_gene_list
        println "INFO: No specific gene list (--gene_list or --imprinted) provided. Plotting will consider all genes from DMRs (if DMRs are found)"
    }

    if (params.plots_only) {
        println "Running in plots-only mode..."
        def annotated_dmrs_input_ch = annotated_dmrs_ch

        boolean use_phased_plotting_logic = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC') && phased_bam_provided
        boolean use_non_phased_plotting_logic = bam_files_provided && !use_phased_plotting_logic

        if (use_phased_plotting_logic) {
            println "INFO: Plotting haplotagged DMRs in plots-only mode."
            if (!params.methylartist_only) {
                plot_phased_dmr_modbamtools(
                    annotated_dmrs_input_ch,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    phased_bam_ch
                )
            }
            if (params.reference) {
                plot_phased_dmr_methylartist(
                    annotated_dmrs_input_ch,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    file(params.reference),
                    phased_bam_ch
                )
            } else if (params.methylartist_only) {
                println "ERROR: --reference genome not provided, but --methylartist_only was specified. Cannot plot haplotagged DMRs with methylartist in plots-only mode."
                exit 1
            } else {
                 println "INFO: --reference not provided, skipping methylartist plots for haplotagged DMRs in plots-only mode."
            }
        } else if (use_non_phased_plotting_logic) {
            println "INFO: Plotting non-phased DMRs in plots-only mode."
            if (!params.methylartist_only) {
                plot_dmr_modbamtools(
                    annotated_dmrs_input_ch,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    input_bam_ch1,
                    input_bam_ch2
                )
            }
            if (params.reference) {
                plot_dmr_methylartist(
                    annotated_dmrs_input_ch,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    file(params.reference),
                    input_bam_ch1,
                    input_bam_ch2
                )
            } else if (params.methylartist_only) {
                println "ERROR: --reference genome not provided, but --methylartist_only was specified. Cannot plot DMRs with methylartist in plots-only mode."
                exit 1
            } else {
                println "INFO: --reference not provided, skipping methylartist plots for DMRs in plots-only mode."
            }
        } else {
            if (phased_bam_provided) {
                 // This case should ideally not be reached if use_phased_plotting_logic was false but phased_bam_provided is true.
                 // It implies flags like --phased_mC were not set.
            } else if (bam_files_provided) {
                 // Similar to above for non-phased.
            }
            // The more general error for missing BAMs is caught by initial validation.
            // If somehow it reaches here, it means logic for selecting plot type failed.
            println "ERROR: In plots-only mode, but suitable BAM files for the specified methylation types are not available or flags are inconsistent."
        }
    } else {
        if (params.'5mC') {
            println "Analysing 5mC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling(group1_beds, group2_beds, 'm', 5, "5mC_group")
            } else {
                prep_bedfile_results1 = prep_bedmethyl_5mC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_5mC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
            }
        } else if (params.'6mA') {
            println "Analysing 6mA methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                dmrs = group_dmr_calling(group1_beds, group2_beds, 'a', 5, "6mA_group")
            } else {  
                prep_bedfile_results1 = prep_bedmethyl_6mA_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_6mA_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)    
            }
        } else if (params.'4mC') {
            println "Analysing 4mC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                // Assuming 4mC uses 'm' and score 5 for group, similar to 5mC group. Adjust if different.
                dmrs = group_dmr_calling(group1_beds, group2_beds, 'm', 5, "4mC_group")
            } else {  
                prep_bedfile_results1 = prep_bedmethyl_4mC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_4mC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)       
            }
        } else if (params.'5hmC') {
            println "INFO: Analysing 5hmC methylation..."
            if (params.input_group1 && params.input_group2) {
                // Group files manually by pattern
                group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
                group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()
                // Assuming 5hmC uses 'h' and score 5 for group. Adjust if different.
                dmrs = group_dmr_calling(group1_beds, group2_beds, 'h', 5, "5hmC_group")
            } else {
                prep_bedfile_results1 = prep_bedmethyl_5hmC_1(input_ch1)
                prep_bedfile_results2 = prep_bedmethyl_5hmC_2(input_ch2)
                dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
            }
        } else if (params.'phased_mC') {
            println "INFO: Analysing haplotagged mC methylation..."
            prep_bedfile_results1 = prep_bedmethyl_mC_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_mC_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else if (params.'phased_hmC') {
            println "INFO: Analysing haplotagged 5hmC methylation..."
            prep_bedfile_results1 = prep_bedmethyl_hmC_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_hmC_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else if (params.'phased_mA') {
            println "INFO: Analysing haplotagged 6mA methylation..."
            prep_bedfile_results1 = prep_bedmethyl_6mA_HP_1(input_ch1)
            prep_bedfile_results2 = prep_bedmethyl_6mA_HP_2(input_ch2)
            dmrs = dmr_calling(prep_bedfile_results1.modified_bed, prep_bedfile_results2.modified_bed)
        } else {
            println "Provide a valid methylation type flag '--5mC, --5hmC, --6mA, --4mC, --phased_mC, --phased_hmC, --phased_mA' and retry, exiting."
            exit 1
        }

        annotated_dmr_output_ch_for_plotting = annotate_dmrs(dmrs.dmr_beds, file(annotationFile)).annotated_dmr_beds
        boolean plot_phased_type = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC')
        boolean is_group_analysis = (params.input_group1 && params.input_group2)
        
        // Only attempt plotting if not a group analysis
        if (is_group_analysis) {
            println "INFO: Plotting is not supported for group analysis."
        } else if (plot_phased_type) {
            if (phased_bam_provided) {
                println "INFO: Plotting haplotagged DMRs."
                if (!params.methylartist_only) {
                    plot_phased_dmr_modbamtools(
                        annotated_dmr_output_ch_for_plotting,
                        gencode_annotation_for_plotting,
                        gene_list_for_plotting,
                        phased_bam_ch
                    )
                }
                if (params.reference) {
                    plot_phased_dmr_methylartist(
                        annotated_dmr_output_ch_for_plotting,
                        gencode_annotation_for_plotting,
                        gene_list_for_plotting,
                        file(params.reference),
                        phased_bam_ch
                    )
                } else if (params.methylartist_only) {
                    println "WARNING: --reference genome not provided, but --methylartist_only was specified. Skipping methylartist plots for haplotagged DMRs."
                } else { // Not methylartist_only and no reference
                    println "INFO: --reference not provided, skipping methylartist plots for haplotagged DMRs."
                }
            } else {
                println "WARNING: Phased BAM file (--phased_modBam) not provided. Skipping plotting for haplotagged DMRs despite phased analysis flag."
            }
        } else if (bam_files_provided) {
            println "INFO: Plotting DMRs."
            if (!params.methylartist_only) {
                plot_dmr_modbamtools(
                    annotated_dmr_output_ch_for_plotting,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    input_bam_ch1,
                    input_bam_ch2
                )
            }
            if (params.reference) {
                plot_dmr_methylartist(
                    annotated_dmr_output_ch_for_plotting,
                    gencode_annotation_for_plotting,
                    gene_list_for_plotting,
                    file(params.reference),
                    input_bam_ch1,
                    input_bam_ch2
                )
            } else if (params.methylartist_only) {
                println "WARNING: --reference genome not provided, but --methylartist_only was specified. Skipping methylartist plots for DMRs."
            } else { // Not methylartist_only and no reference
                 println "INFO: --reference not provided, skipping methylartist plots for DMRs."
            }
        } else {
            println "INFO: No suitable BAM files provided (standard or phased). Skipping plotting of DMRs."
        }    
    }
}