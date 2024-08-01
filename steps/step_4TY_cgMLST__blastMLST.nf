nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;getGenusSpeciesOptionalUnkeyed } from '../functions/parameters.nf'
include { stepInputs;isSpeciesSupported } from '../functions/common.nf'

def schemaBlastMLST = "cgMLST_lmono_Pasteur_20191007"
def schemaBlastMLSTPath = "/mnt/biowork/databases/PROGRAMS/blastMLST/cgMLST_lmono_Pasteur_20191007.tgz"

def lociSchemaPasteur = "/mnt/biowork/databases/PROGRAMS/chewbbaca/L_monocytogenes_schemas/Pasteur_1748/cgMLST_20191007.txt"

def ex = executionMetadata()

def STEP = '4TY_cgMLST'
def METHOD = 'blastMLST' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_SPECIES_ALLOWED = [
  "Listeria_monocytogenes"
]

process blast_cgmlst {
    container "${LOCAL_REGISTRY}/bioinfo/mlst:2.16.1--7b565d3e3f"
    containerOptions = "-v /mnt/biowork:/mnt/biowork:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro  -u 0:0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    when:
      isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, scaffolds200, task.process)        
    input:
      tuple val(riscd_input), path(scaffolds200)
      val genus_species
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_blast_cgmlst.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        /scripts/build-db.sh ${schemaBlastMLSTPath}
        mlst --threads 16 --scheme ${schemaBlastMLST} ${scaffolds200} &> ${base}_blastMLST.log
        /scripts/process-blast-result.py ${md.ds} ${base}_blastMLST.log ${lociSchemaPasteur} ${base}_blastMLST_results_alleles.tsv
      """
}

workflow step_4TY_cgMLST__blastMLST {
    take: 
      assembly
      speciesCode
    main:
      //assume channels are already crossed
      blast_cgmlst(assembly,speciesCode)
}

workflow {
    step_4TY_cgMLST__blastMLST(getInput(), getGenusSpeciesOptionalUnkeyed())
}
