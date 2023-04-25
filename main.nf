#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

// Sample sheet (TSV file)
// Columns: name, highlight, txt
//   name: Name of the sample to be applied as the title of graph (string)
//   highlight: Which, if any, chromosomal arms to highlight (CSV list)
//   txt: SMASH results file (TSV file) -- these files are looked for in `data/`
//
// Example:
// name    highlight       txt
// A2780 TK+ AAVS1 gRNA c1 1q      C2-Rosa-03_20000_bins_results.txt
// A2780 TK+ TP53 gRNA c1  1q      C2-286-04_20000_bins_results.txt
// A2780 TK+ TP53 gRNA c2  1q      C2-287-05_20000_bins_results.txt
// ...
params.sampleSheet = "SMASH-Graphs.txt"

// Definition of chromosomal arms breaks/lengths (CSV file)
params.arms = "${projectDir}/ref/hg19-chromosome-arm-lengths.csv"

// The vizualization extensions to produce (comma-separated list)
params.exts = "png,pdf"

workflow {
    // Get gene set from file into Groovy/NF List for reuse where relevant
    def samples = Channel.fromPath(params.sampleSheet, checkIfExists: true) \
                | splitCsv(header: true, sep: '\t') \
                | map { row -> tuple(row."name", row."highlight", file("data/${row.txt}") ) }
    def arms = Channel.fromPath(params.arms, checkIfExists: true).first()
    def exts = Channel.fromList(params.exts?.tokenize(','))

    produce_viz(samples, arms, exts)
}

process produce_viz {
    conda 'conda-forge::r-tidyverse=1.3.2 conda-forge::r-argparser=0.7.1 conda-forge::r-stringr=1.4.1'
    label 'process_medium'
    tag "${cleanname}"
    publishDir "results/${ext}/", mode: 'copy'

    input:
      tuple val(name),
            val(highlight),
            path(txt)
      path arms
      each ext
    
    output:
      path "${cleanname}.${ext}"

    script:
      cleanname = name.replaceAll('[^a-zA-Z0-9_]+', '_')
      def optional = highlight != '' ? "--highlight ${highlight}" : ''
      """
      smash-viz.R \
        --arms ${arms} \
        --input ${txt} \
        ${optional} \
        --color blue \
        --title '${name}' \
        --output '${cleanname}.${ext}'
      """
}
