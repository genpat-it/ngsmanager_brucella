nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '0SQ_rawreads'
def METHOD = 'tenderECDC' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process tender_ecdc_qc {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3500.MB, task.attempt ) }
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_input2), path(genusReport)
      tuple val(_), val(species)
    output:
      path '{*.tsv,*.json}'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_input2], md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'fail_QC.tsv', saveAs: { "${base}_fail_QC.tsv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_tenderECDC_md5qc.tsv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_preprocessing.cfg" }
    script:
      (r1,r2) = reads
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      /scripts/tenderQC_report.py -n ${md.cmp} -r1 ${r1} -r2 ${r2} -g ${genusReport} -sp '${species}' > ${base}_tenderECDC_md5qc.tsv
      """
}

workflow step_0SQ_rawreads__tenderECDC {
  take: 
      reads
      genus_report
      species    
  main:
      tender_ecdc_qc(reads, genus_report, species)
}