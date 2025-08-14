#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Help message
def helpMessage = """
Usage: nextflow run main.nf [options]

Options:
  --input_file1        First input bedmethyl file from sample_1 or haplotype_1 bedmethyl (required for analysis)
  --input_file2        Second input bedmethyl file from sample_2 or haplotype_2 bedmethyl (required for analysis)
  --output_dir         Output directory (required)
  --input_modBam1      First input modified BAM used to generate bedmethyl for sample_1 (optional, for plotting or pileup)
  --input_modBam2      Second input modified BAM used to generate bedmethyl for sample_2 (optional, for plotting or pileup)
  --input_group1       Directory containing bedmethyl files (with *.bed extension) for group 1 (optional, for group analysis)
  --input_group2       Directory containing bedmethyl files (with *.bed extension) for group 2 (optional, for group analysis)
  --gene_list          A list of genes to generate plots on (optional). If not provided and --imprinted is used, imprinted genes will be plotted.
  --reference          Reference genome FASTA file (required for methylartist plots and pileup)
  --imprinted          Use this flag to plot DMRs overlapping imprinted genes if --gene_list is not provided (optional, for haplotagged plots)
  --5mC                Use this flag to trigger 5mC DMR calling (optional)
  --5hmC               Use this flag to trigger 5hmC DMR calling (optional)
  --6mA                Use this flag to trigger 6mA DMR calling (optional)
  --4mC                Use this flag to trigger 4mC DMR calling (optional)
  --phased_mC          Use this flag to trigger haplotagged 5mC/4mC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_mA          Use this flag to trigger haplotagged 6mA DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_hmC         Use this flag to trigger haplotagged 5hmC DMR calling if input files are haplotagged bedmethyls (optional)
  --phased_modBam      Haplotagged modified BAM (required for plotting DMRs if --phased_mC, --phased_mA or --phased_hmC is used)
  --pileup             Use this flag to call methylated bases from modified BAMs using modkit (requires --reference and relevant BAM files)
  --plot               Enable plotting after DMR calling (off by default; works with pileup or bedmethyl inputs)
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
params.pileup = false
params.plot = false

// Help information
if (params.help) {
    println helpMessage
    exit 0
}

// Validation for pileup mode
if (params.pileup) {
    if (!params.reference) {
        println "ERROR: --reference genome is required when --pileup is specified."
        println helpMessage
        exit 1
    }
    
    // Check for phased pileup
    boolean phased_pileup = (params.'phased_mC' || params.'phased_mA') && params.phased_modBam
    boolean standard_pileup = (params.'5mC' || params.'6mA') && params.input_modBam1 && params.input_modBam2
    
    if (!phased_pileup && !standard_pileup) {
        println "ERROR: For --pileup mode, provide either:"
        println "  - Standard pileup: --input_modBam1, --input_modBam2 with --5mC or --6mA"
        println "  - Phased pileup: --phased_modBam with --phased_mC or --phased_mA"
        println helpMessage
        exit 1
    }
}

// Early validation for --reference if --methylartist_only is used (only when plotting)
if ((params.plot || params.plots_only) && params.methylartist_only && !params.reference) {
    println "ERROR: --reference genome is required when --methylartist_only is specified for plotting."
    println helpMessage
    exit 1
}

// Separate validation for plots_only mode
if (params.plots_only) {
    if (!params.annotated_dmrs || !params.output_dir) {
        println "ERROR: In plots-only mode, --annotated_dmrs and --output_dir are required."
        println helpMessage
        exit 1
    }
    boolean phased_bams_for_plots = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC') && params.phased_modBam
    boolean non_phased_bams_for_plots = params.input_modBam1 && params.input_modBam2
    if (!phased_bams_for_plots && !non_phased_bams_for_plots) {
        println "ERROR: In plots-only mode, relevant BAM files (--phased_modBam or --input_modBam1/2) are required for plotting."
        println helpMessage
        exit 1
    }
} else if (!params.pileup) {
    // Validation for regular analysis mode (not pileup, not plots_only)
    if ((!params.input_file1 && !params.input_group1) || 
        (!params.input_file2 && !params.input_group2) || 
        !params.output_dir) {
        println "ERROR: In analysis mode, input files (--input_file1/2 or --input_group1/2) and --output_dir are required"
        println helpMessage
        exit 1
    }
}

// Check if output directory exists and create one if not already present
def outputDir = file(params.output_dir)
if (!outputDir.exists()) {
    println "INFO: Output directory does not exist. Creating one: ${params.output_dir}"
    outputDir.mkdirs()
}

// Input channels
input_ch1 = params.input_file1 ? Channel.fromPath(params.input_file1, checkIfExists: true) : Channel.empty()
input_ch2 = params.input_file2 ? Channel.fromPath(params.input_file2, checkIfExists: true) : Channel.empty()
annotationFile = file("${workflow.projectDir}/annotations/gencode.v46.annotation.exon-promoters-introns.sorted.bed") // For annotate_dmrs

// Input channel for pileup reference
pileup_ref_ch = params.pileup ? Channel.fromPath(params.reference, checkIfExists: true).map { [it.baseName, file(it), file(it + '.fai')] } : Channel.empty()

// Reference channel for plotting (independent of pileup)
plot_ref_ch = params.reference ? Channel.fromPath(params.reference, checkIfExists: true).map { [it.baseName, file(it), file(it + '.fai')] } : Channel.empty()

// Ensure FASTA index exists if plotting needs methylartist
if ((params.plot || params.plots_only || params.methylartist_only) && params.reference && !file("${params.reference}.fai").exists()) {
    println "ERROR: Reference index not found: ${params.reference}.fai. Run: samtools faidx ${params.reference}"
    exit 1
}

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
if (!params.gene_list && !params.imprinted) {
    empty_gene_list.text = ""
}

// Include processes
include { prep_bedmethyl_5mC  as prep_bedmethyl_5mC_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5mC  as prep_bedmethyl_5mC_2  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5hmC as prep_bedmethyl_5hmC_1 } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5hmC as prep_bedmethyl_5hmC_2 } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_6mA  as prep_bedmethyl_6mA_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_6mA  as prep_bedmethyl_6mA_2  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5mC  as prep_bedmethyl_4mC_1  } from './workflows/prep_bedmethyl.nf'
include { prep_bedmethyl_5mC  as prep_bedmethyl_4mC_2  } from './workflows/prep_bedmethyl.nf'
include { prep_phased_bedmethyl_mC_HP  as prep_bedmethyl_mC_HP_1  } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mC_HP  as prep_bedmethyl_mC_HP_2  } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_hmC_HP as prep_bedmethyl_hmC_HP_1 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_hmC_HP as prep_bedmethyl_hmC_HP_2 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mA_HP  as prep_bedmethyl_6mA_HP_1 } from './workflows/prep_phased_bedmethyl.nf'
include { prep_phased_bedmethyl_mA_HP  as prep_bedmethyl_6mA_HP_2 } from './workflows/prep_phased_bedmethyl.nf'
include { modkit_phased_mC; modkit_phased_mA } from './workflows/modkit.nf'
include { modkit_mC as modkit_mC_1 } from './workflows/modkit.nf'
include { modkit_mC as modkit_mC_2 } from './workflows/modkit.nf'
include { modkit_mA as modkit_mA_1 } from './workflows/modkit.nf'
include { modkit_mA as modkit_mA_2 } from './workflows/modkit.nf'
include { dmr_calling } from './workflows/dmr_calling.nf'
include { split_group_beds_by_chr as split_beds_group1 } from './workflows/dmr_calling.nf'
include { split_group_beds_by_chr as split_beds_group2 } from './workflows/dmr_calling.nf'
include { group_dmr_calling } from './workflows/dmr_calling.nf'
include { aggregate_dmrs } from './workflows/aggregate_dmrs.nf'
include { annotate_dmrs } from './workflows/annotate_dmrs.nf'
include { report_dmrs } from './workflows/generate_dmr_report.nf'
include { plot_dmr_modbamtools; plot_dmr_methylartist; plot_phased_dmr_modbamtools; plot_phased_dmr_methylartist } from './workflows/plot_dmrs.nf'

workflow {
    // Determine the gene list to use for plotting
    def gene_list_for_plotting = getGeneListForPlotting()

    if (params.plots_only) {
        runPlotsOnly(gene_list_for_plotting)
    } else if (params.pileup) {
        runPileupMode(gene_list_for_plotting)
    } else {
        runAnalysisMode(gene_list_for_plotting)
    }
}

def getGeneListForPlotting() {
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
                println "INFO: Using processed imprinted gene list for plotting or annotation: ${processed_file.name}"
                return processed_file
            } catch (Exception e) {
                println "ERROR: Failed to process imprinted genes file ${imprinted_file.name}: ${e.getMessage()}"
                exit 1
            }
        } else {
            println "ERROR: --imprinted flag was used, but the imprinted genes file was not found."
            exit 1
        }
    } else if (params.gene_list) {
        println "INFO: Using user-provided gene list for plotting: ${params.gene_list}"
        return file(params.gene_list)
    } else {
        if (!empty_gene_list.exists()) {
            empty_gene_list.text = ""
        }
        println "INFO: No specific gene list (--gene_list or --imprinted) provided. Plotting will consider all genes from DMRs (if DMRs are found)"
        return empty_gene_list
    }
}

def runPlotsOnly(gene_list_for_plotting) {
    println "MODE: Plots-only..."
    def annotated_dmrs_input_ch = annotated_dmrs_ch

    boolean use_phased_plotting_logic = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC') && phased_bam_provided
    boolean use_non_phased_plotting_logic = bam_files_provided && !use_phased_plotting_logic

    if (use_phased_plotting_logic) {
        println "INFO: Plotting haplotagged DMRs in plots-only mode."
        plotPhasedDMRs(annotated_dmrs_input_ch, gene_list_for_plotting, phased_bam_ch)
    } else if (use_non_phased_plotting_logic) {
        println "INFO: Plotting non-phased DMRs in plots-only mode."
        plotStandardDMRs(annotated_dmrs_input_ch, gene_list_for_plotting, input_bam_ch1, input_bam_ch2)
    } else {
        println "ERROR: In plots-only mode, but suitable BAM files for the specified methylation types are not available or flags are inconsistent."
    }
}

def runPileupMode(gene_list_for_plotting) {
    println "MODE: Pileup..."
    
    if (params.'phased_mC' && params.phased_modBam) {
        println "INFO: Running phased 5mC pileup and analysis..."
        phased_mC_beds = modkit_phased_mC(phased_bam_ch, pileup_ref_ch)
        
        // Process the phased beds through prep_phased_bedmethyl
        hp1_beds_ch = phased_mC_beds.modkit_H1.map { name, bed -> bed }
        hp2_beds_ch = phased_mC_beds.modkit_H2.map { name, bed -> bed }
        
        prep_ch1 = prep_bedmethyl_mC_HP_1(hp1_beds_ch)
        prep_ch2 = prep_bedmethyl_mC_HP_2(hp2_beds_ch)
        
        dmrs = runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
        annotated_dmr_output = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
        report_dmrs(annotated_dmr_output.annotated_dmr_beds, annotated_dmr_output.annotation_log)
        
        // Plot if requested and BAM provided
        if (params.plot && phased_bam_provided) {
            plotPhasedDMRs(annotated_dmr_output.annotated_dmr_beds, gene_list_for_plotting, phased_bam_ch)
        } else if (!params.plot) {
            println "INFO: --plot not specified; skipping plotting."
        }
        
    } else if (params.'phased_mA' && params.phased_modBam) {
        println "INFO: Running phased 6mA pileup and analysis..."
        phased_mA_beds = modkit_phased_mA(phased_bam_ch, pileup_ref_ch)
        
        // Process the phased beds through prep_phased_bedmethyl
        hp1_beds_ch = phased_mA_beds.modkit_H1.map { name, bed -> bed }
        hp2_beds_ch = phased_mA_beds.modkit_H2.map { name, bed -> bed }
        
        prep_ch1 = prep_bedmethyl_6mA_HP_1(hp1_beds_ch)
        prep_ch2 = prep_bedmethyl_6mA_HP_2(hp2_beds_ch)
        
        dmrs = runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
        annotated_dmr_output = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
        report_dmrs(annotated_dmr_output.annotated_dmr_beds, annotated_dmr_output.annotation_log)
        
        // Plot if requested and BAM provided
        if (params.plot && phased_bam_provided) {
            plotPhasedDMRs(annotated_dmr_output.annotated_dmr_beds, gene_list_for_plotting, phased_bam_ch)
        } else if (!params.plot) {
            println "INFO: --plot not specified; skipping plotting."
        }
        
    } else if (params.'5mC' && params.input_modBam1 && params.input_modBam2) {
        println "INFO: Running standard 5mC pileup and analysis..."
        mC_bed1 = modkit_mC_1(input_bam_ch1, pileup_ref_ch)
        mC_bed2 = modkit_mC_2(input_bam_ch2, pileup_ref_ch)
        
        // Process the beds through prep_bedmethyl
        bed1_ch = mC_bed1.modkit.map { name, bed -> bed }
        bed2_ch = mC_bed2.modkit.map { name, bed -> bed }
        
        prep_ch1 = prep_bedmethyl_5mC_1(bed1_ch)
        prep_ch2 = prep_bedmethyl_5mC_2(bed2_ch)
        
        dmrs = runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
        annotated_dmr_output = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
        report_dmrs(annotated_dmr_output.annotated_dmr_beds, annotated_dmr_output.annotation_log)
        
        // Plot if requested and BAMs provided
        if (params.plot && bam_files_provided) {
            plotStandardDMRs(annotated_dmr_output.annotated_dmr_beds, gene_list_for_plotting, input_bam_ch1, input_bam_ch2)
        } else if (!params.plot) {
            println "INFO: --plot not specified; skipping plotting."
        }
        
    } else if (params.'6mA' && params.input_modBam1 && params.input_modBam2) {
        println "INFO: Running standard 6mA pileup and analysis..."
        mA_bed1 = modkit_mA_1(input_bam_ch1, pileup_ref_ch)
        mA_bed2 = modkit_mA_2(input_bam_ch2, pileup_ref_ch)
        
        // Process the beds through prep_bedmethyl
        bed1_ch = mA_bed1.modkit.map { name, bed -> bed }
        bed2_ch = mA_bed2.modkit.map { name, bed -> bed }
        
        prep_ch1 = prep_bedmethyl_6mA_1(bed1_ch)
        prep_ch2 = prep_bedmethyl_6mA_2(bed2_ch)
        
        dmrs = runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
        annotated_dmr_output = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
        report_dmrs(annotated_dmr_output.annotated_dmr_beds, annotated_dmr_output.annotation_log)
        
        // Plot if requested and BAMs provided
        if (params.plot && bam_files_provided) {
            plotStandardDMRs(annotated_dmr_output.annotated_dmr_beds, gene_list_for_plotting, input_bam_ch1, input_bam_ch2)
        } else if (!params.plot) {
            println "INFO: --plot not specified; skipping plotting."
        }
    }
}

def runAnalysisMode(gene_list_for_plotting) {
    // Validate methylation type flags
    def methylation_types = [params.'5mC', params.'5hmC', params.'6mA', params.'4mC', 
                           params.'phased_mC', params.'phased_hmC', params.'phased_mA']
    if (!methylation_types.any()) {
        println "ERROR: Provide a valid methylation type flag '--5mC, --5hmC, --6mA, --4mC, --phased_mC, --phased_hmC, --phased_mA' and retry."
        exit 1
    }

    // Run analysis based on methylation type
    def dmrs
    if (params.'5mC') {
        println "Analysing 5mC methylation..."
        dmrs = runMethylationAnalysis('5mC', { prep_bedmethyl_5mC_1(it) }, { prep_bedmethyl_5mC_2(it) })
    } else if (params.'6mA') {
        println "Analysing 6mA methylation..."
        dmrs = runMethylationAnalysis('6mA', { prep_bedmethyl_6mA_1(it) }, { prep_bedmethyl_6mA_2(it) })
    } else if (params.'4mC') {
        println "Analysing 4mC methylation..."
        dmrs = runMethylationAnalysis('4mC', { prep_bedmethyl_4mC_1(it) }, { prep_bedmethyl_4mC_2(it) })
    } else if (params.'5hmC') {
        println "INFO: Analysing 5hmC methylation..."
        dmrs = runMethylationAnalysis('5hmC', { prep_bedmethyl_5hmC_1(it) }, { prep_bedmethyl_5hmC_2(it) })
    } else if (params.'phased_mC') {
        println "INFO: Analysing haplotagged mC methylation..."
        dmrs = runPhasedAnalysis({ prep_bedmethyl_mC_HP_1(it) }, { prep_bedmethyl_mC_HP_2(it) })
    } else if (params.'phased_hmC') {
        println "INFO: Analysing haplotagged 5hmC methylation..."
        dmrs = runPhasedAnalysis({ prep_bedmethyl_hmC_HP_1(it) }, { prep_bedmethyl_hmC_HP_2(it) })
    } else if (params.'phased_mA') {
        println "INFO: Analysing haplotagged 6mA methylation..."
        dmrs = runPhasedAnalysis({ prep_bedmethyl_6mA_HP_1(it) }, { prep_bedmethyl_6mA_HP_2(it) })
    }

    // Annotate and report DMRs
    annotated_dmr_output = annotate_dmrs(dmrs.dmr_beds, file(annotationFile))
    report_dmrs(annotated_dmr_output.annotated_dmr_beds, annotated_dmr_output.annotation_log)
    
        // Plot DMRs if not group analysis

    boolean is_group_analysis = (params.input_group1 && params.input_group2)
    if (!is_group_analysis) {
        if (params.plot) {
            plotDMRsBasedOnType(annotated_dmr_output.annotated_dmr_beds, gene_list_for_plotting)
        } else {
            println "INFO: --plot not specified; skipping plotting."
        }
    } else {
        if (params.plot) {
            println "INFO: Plotting is not supported for group analysis."
        }
    }    
}

def runMethylationAnalysis(methylation_type, prep_func1, prep_func2) {
    if (params.input_group1 && params.input_group2) {
        return runGroupAnalysis(methylation_type)
    } else {
        return runStandardAnalysis(prep_func1, prep_func2)
    }
}

def runGroupAnalysis(methylation_type) {
    group1_beds = Channel.fromPath("${params.input_group1}/*.bed").collect()
    group2_beds = Channel.fromPath("${params.input_group2}/*.bed").collect()

    // Determine parameters based on methylation type
    def base_char = (methylation_type == '6mA') ? 'a' : 
                   (methylation_type == '5hmC') ? 'h' : 'm'
    def position = 5

    group1_split = split_beds_group1(group1_beds, "group1", base_char, position)
    group2_split = split_beds_group2(group2_beds, "group2", base_char, position)

    group1_by_chr = group1_split.chr_beds.flatten()
        .map { file -> [file.name.split('_')[0], file] }
        .groupTuple()

    group2_by_chr = group2_split.chr_beds.flatten()
        .map { file -> [file.name.split('_')[0], file] }
        .groupTuple()

    chr_grouped = group1_by_chr.join(group2_by_chr, by: 0)
    chr_dmrs = group_dmr_calling(chr_grouped, "${methylation_type}_group")

    all_dmr_files = chr_dmrs.chr_dmrs.collect()
    all_status_logs = chr_dmrs.status_log.collect()
    all_debug_dirs = chr_dmrs.debug_output.collect()

    return aggregate_dmrs(all_dmr_files, all_status_logs, all_debug_dirs)
}

def runStandardAnalysis(prep_func1, prep_func2) {
    prep_ch1 = prep_func1(input_ch1)
    prep_ch2 = prep_func2(input_ch2)
    return runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
}

def runPhasedAnalysis(prep_func1, prep_func2) {
    prep_ch1 = prep_func1(input_ch1)
    prep_ch2 = prep_func2(input_ch2)
    return runDMRAnalysisFromPrepped(prep_ch1, prep_ch2)
}

def runDMRAnalysisFromPrepped(prep_ch1, prep_ch2) {
    chr_files_ch1 = prep_ch1.chr_beds.flatten()
        .map { file -> [file.baseName.split('_')[0], file] }

    chr_files_ch2 = prep_ch2.chr_beds.flatten()
        .map { file -> [file.baseName.split('_')[0], file] }

    chr_pairs = chr_files_ch1.join(chr_files_ch2, by: 0)
    chr_dmrs = dmr_calling(chr_pairs)

    all_dmr_files = chr_dmrs.chr_dmrs.collect()
    all_status_logs = chr_dmrs.status_log.collect()
    all_debug_dirs = chr_dmrs.debug_output.collect()

    return aggregate_dmrs(all_dmr_files, all_status_logs, all_debug_dirs)
}

def plotDMRsBasedOnType(annotated_dmr_output_ch, gene_list_for_plotting) {
    boolean plot_phased_type = (params.'phased_mC' || params.'phased_mA' || params.'phased_hmC')

    if (plot_phased_type && phased_bam_provided) {
        plotPhasedDMRs(annotated_dmr_output_ch, gene_list_for_plotting, phased_bam_ch)
    } else if (bam_files_provided) {
        plotStandardDMRs(annotated_dmr_output_ch, gene_list_for_plotting, input_bam_ch1, input_bam_ch2)
    } else {
        println "INFO: No suitable BAM files provided for plotting."
    }
}

def plotPhasedDMRs(dmr_ch, gene_list, phased_bam_ch) {
    println "INFO: Plotting haplotagged DMRs..."
    if (!params.methylartist_only) {
        plot_phased_dmr_modbamtools(dmr_ch, gencode_annotation_for_plotting, gene_list, phased_bam_ch)
    }
    if (params.reference) {
        plot_phased_dmr_methylartist(dmr_ch, gencode_annotation_for_plotting, gene_list, plot_ref_ch, phased_bam_ch)
    } else if (params.methylartist_only) {
        println "WARNING: --reference genome not provided, but --methylartist_only was specified. Skipping methylartist plots for haplotagged DMRs."
    } else {
        println "INFO: --reference not provided, skipping methylartist plots for haplotagged DMRs."
    }
}

def plotStandardDMRs(dmr_ch, gene_list, bam_ch1, bam_ch2) {
    println "INFO: Plotting DMRs..."
    if (!params.methylartist_only) {
        plot_dmr_modbamtools(dmr_ch, gencode_annotation_for_plotting, gene_list, bam_ch1, bam_ch2)
    }
    if (params.reference) {
        plot_dmr_methylartist(dmr_ch, gencode_annotation_for_plotting, gene_list, plot_ref_ch, bam_ch1, bam_ch2)
    } else if (params.methylartist_only) {
        println "WARNING: --reference genome not provided, but --methylartist_only was specified. Skipping methylartist plots for DMRs."
    } else {
        println "INFO: --reference not provided, skipping methylartist plots for DMRs."
    }
}