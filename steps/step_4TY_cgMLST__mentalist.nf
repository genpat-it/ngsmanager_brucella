nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getSingleInput;getGenusSpeciesOptionalUnkeyed;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def dbMentalistListeria="/mnt/biowork/databases/PROGRAMS/mentalist/l_mono_Pasteur_cgMLST_20191007_k21.db"

def ex = executionMetadata()

def STEP = '4TY_cgMLST'
def METHOD = 'mentalist' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_SPECIES_ALLOWED = [
  "Listeria_monocytogenes"
]

process mentalist {
    container "${LOCAL_REGISTRY}/bioinfo/mentalist:1.0.0--39e9e05e54"
    containerOptions = "-v /mnt/biowork:/mnt/biowork:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 5.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, reads, task.process)  && isCompatibleWithSeqType(reads, ['illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      val genus_species
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_mentalist*.*'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_mentalist.cfg" }
    script:
      (t1,t2,u) = reads
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        mentalist call -o ./${base}_mentalist --output_votes --output_special --db ${dbMentalistListeria} -1 ${t1} -2 ${t2} 2> ${base}_mentalist.log 
        /scripts/process-mentalist-result.py ${md.ds} ${base}_mentalist ${base}_mentalist_results_alleles.tsv
      """
}

workflow step_4TY_cgMLST__mentalist {
    take: 
      trimmed
      speciesCode
    main:
      //assume channels are already crossed
      mentalist(trimmed,speciesCode)
}

workflow {
    step_4TY_cgMLST__mentalist(getSingleInput(), getGenusSpeciesOptionalUnkeyed())
}