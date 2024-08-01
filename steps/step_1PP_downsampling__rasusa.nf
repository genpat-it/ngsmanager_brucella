nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getSingleInput;optionalOption;optionalOptionWithKey;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_downsampling'
def METHOD = 'rasusa' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"
def PARAMS_PREFIX = "${STEP}__${METHOD}___"

process rasusa {
    container "quay.io/biocontainers/rasusa:2.0.0--h031d066_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"   
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [opts:opts])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*_input.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'rasusa.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_rasusa"

      opts = optionalOption(PARAMS_PREFIX, 'num')
      opts += optionalOption(PARAMS_PREFIX, 'bases')
      opts += optionalOptionWithKey(PARAMS_PREFIX, 'genome-size', 'genome_size')
      opts += optionalOption(PARAMS_PREFIX, 'coverage')
      opts += optionalOption(PARAMS_PREFIX, 'frac') 

       if (r2) { 
        """
        rasusa reads ${opts} -s 1 ${r1} ${r2} -o ${base}_R1.fastq.gz -o ${base}_R2.fastq.gz -v 2> rasusa.log
        """
      } else {
        """
        rasusa reads ${opts} -s 1 ${r1} -o ${base}.fastq.gz -v 2> rasusa.log
        """        
      }     
}

workflow step_1PP_downsampling__rasusa {
    take: 
      reads    
    main:
      rasusa(reads)
}

workflow {
    step_1PP_downsampling__rasusa(getSingleInput())
}