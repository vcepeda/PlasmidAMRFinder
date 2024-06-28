#!/usr/bin/env nextflow

/* 
===============================================================================
                              P L A S M I D A M R F I N D E R  
===============================================================================
Nextflow pipeline for resistance plasmid identification and annotation using 
Nanopore reads and bacterial genome assemblies
-------------------------------------------------------------------------------
@ Author
Victoria Cepeda 
Caspar Groß (adapted from PLASMIDIDENT)
-------------------------------------------------------------------------------
@ Documentation
https://github.com/vcepeda/PlasmidAMRFinder/README.md
------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

workflow {
    plasmidAmrFinder()
}

workflow plasmidAmrFinder {
    // Check special input parameters
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.version) {
        pipelineMessage()
        exit 0
    }
    if (!params.input) {
        helpMessage()
        exit 0
    }

    // Setup
    samples = getFiles(params.input)
    env = params.env
    startMessage()
    runParamCheck()

    // Execute the processes
    samples | filter_reads | view
}

process filter_reads {
    tag { id }

    input:
    tuple val(id), path(assembly), path(lr) from samples

    output:
    tuple val(id), path(assembly), path('reads_filtered.fastq')

    script:
    if (params.noSubsampling) {
        """
        zcat -f ${lr} > reads_filtered.fastq
        """
    } else {
        """
        ${env}
        len=\$(grep -v '>' ${assembly} | wc -c)
        nbases=\$(expr \$len * ${params.mappingCov})
        filtlong -t \$nbases --length_weight 0 ${lr} > reads_filtered.fastq
        """
    }
}

def asciiArt() {
    log.info "         _                     _     _      _    __  __ ____     __ _           _           "
    log.info "  _ __ | | __ _ ___ _ __ ___ (_) __| |    / \\  |  \\/  |  _ \\   / _(_)_ __   __| | ___ _ __ "
    log.info " | '_ \\| |/ _` / __| '_ ` _ \\| |/ _` |   / _ \\ | |\\/| | |_) | | |_| | '_ \\ / _` |/ _ \\ '__|"
    log.info " | |_) | | (_| \\__ \\ | | | | | | (_| |  / ___ \\| |  | |  _ <  |  _| | | | | (_| |  __/ |   "
    log.info " | .__/|_|\\__,_|___/_| |_| |_|_|\\__,_| /_/   \\_\\_|  |_|_| \\_\\ |_| |_|_| |_|\\__,_|\\___|_|   "
    log.info " |_|                                                                                       "
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def runParamCheck() {
    if (params.minLength < params.seqPadding) {
        log.info "Minimum contig length cannot be shorter than seqPadding. Please adjust parameters and restart"
    }
}

def getFiles(tsvFile) {
    // Extracts Read Files from TSV
    if (workflow.profile.contains('test')) {
        inputFile = file("$baseDir/" + tsvFile)
    } else {
        inputFile = file(tsvFile)
    }
    log.info "------------------------------"
    Channel.fromPath(inputFile)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .splitCsv(sep: '\t')
        .map { row ->
            [row[0], file(row[1]), file(row[2])]
        }
}

def helpMessage() {
    // Display help message
    log.info """  Usage:
       nextflow run main.nf --input <file.csv> [options]
    --input <file.tsv>
       TSV file containing paths to files (id | assembly | longread)
  Parameters:
    --outDir
    Output location (Default: current working directory)
    --maxLength <bases> (Default: 500000)
    Contigs larger than maxLength will not be considered a putative plasmid
    --seqPadding <bases> (Default: 2000)
    Length of recycled sequences at contig edges for long read mapping.
    --covWindow <bases> (Default: 50)
    Moving window size for coverage calculation
    --mappingCov <coverage> (Default: 50)
    Target coverage for long read sampling
    --noSubsampling
    Skips the read subsampling step. Use when read coverage is not uniform.
    --version
      Displays pipeline version
    --help
      Displays this help
  Profiles:
    -profile local
    Pipeline runs with locally installed conda environments (found in env/ folder)
    -profile test
    Runs complete pipeline on small included test dataset
    -profile localtest
    Runs test profile with locally installed conda environments
    """
}

def grabRevision() {
    // Return the same string executed from github or not
    return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage() {
    // Minimal information message
    log.info "Command Line  : " + workflow.commandLine
    log.info "Input file    : " + params.inputFile
    log.info "Profile       : " + workflow.profile
    log.info "Project Dir   : " + workflow.projectDir
    log.info "Launch Dir    : " + workflow.launchDir
    log.info "Work Dir      : " + workflow.workDir
    log.info "Cont Engine   : " + workflow.containerEngine
    log.info "Out Dir       : " + params.outDir
    log.info "Align. Overlp.: " + params.seqPadding
    log.info "Cov. window   : " + params.covWindow
    log.info "Max Plasm. Len: " + params.maxLength
    log.info "Min Plasm. Len: " + params.minLength
    log.info "Target cov.   : " + params.mappingCov
    log.info "read sampling : " + !params.noSubsampling
    log.info "Containers    : " + workflow.container 
}

def pipelineMessage() {
    // Display hybridAssembly info message
    log.info "PlasmidAMRFinder Pipeline ~  version ${workflow.manifest.version} - revision " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
    // Display start message
    this.asciiArt()
    this.pipelineMessage()
    this.minimalInformationMessage()
}

workflow.onComplete {
    // Display complete message
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
}
