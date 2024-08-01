nextflow.enable.dsl=2

include { getInput;param } from '../functions/parameters.nf'
include { flattenPath } from '../functions/common.nf'

def REFERENCE_PATH = "${params.assets_dir}/multi_alignment__gisaid/NC_045512.fasta"

process gisaid_preprocess {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/multi_alignment__gisaid:/scripts:ro"
    input:
      path(fastas)
      path(aliases)
      path(reference)
    output:
      path('consensus.fasta'), emit : multifasta
      path("*")
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log', saveAs: { "gisaid_preprocess.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "gisaid_preprocess.cfg" }
    script:
      """
      #!/bin/bash -euo pipefail
      ls ${fastas} >> fasta_list.txt
      /scripts/prepare-gisaid-input.py --file_list fasta_list.txt --aliases ${aliases} --reference ${reference} --output consensus.fasta
      """
}


process muscle {
    container "quay.io/biocontainers/muscle:3.8.1551--h7d875b9_6"
    input:
      path(multifasta)
    output:
      path("*")
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log', saveAs: { "muscle.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "muscle.cfg" }
    script:
      """
        muscle -in ${multifasta} -out allineamento.fasta -maxiters 1 -diags
      """
}

workflow multi_alignment__gisaid {
    take: 
        fastas
        aliases
    main:      
        gisaid_preprocess(fastas, aliases, file(REFERENCE_PATH)).multifasta | muscle
}

workflow {
    getInput()
        .map { it[1] }
        .collect()
        .set { fastas }      
    multi_alignment__gisaid(fastas, param('aliases'));
}
