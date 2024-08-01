nextflow.enable.dsl=2

include { extractKey;taskMemory;executionMetadata;parseMetadataFromFileName } from '../functions/common.nf'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { getSingleInput;getReferences } from '../functions/parameters.nf'

process leggitarget_process {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/module_LEO:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      tuple path(samtools_depth), val(reference)
    output:
      path '*.sh', hidden: true
      tuple path("${md.cmp}_LEO_match.txt"), val(reference), emit: target_match   
      publishDir mode: 'rellink', "${params.outdir}", pattern: '*.txt'
    script:
      md = parseMetadataFromFileName(samtools_depth.getName())
      """
        /scripts/leggimatch.py ${md.cmp} ${md.ds} ${samtools_depth} ${md.cmp}_LEO_match.txt
	    """
}

workflow module_LEO {
    take: 
        reads
        reference 
    main:
        coverage_depth  = step_2AS_mapping__ivar(reads, reference).coverage_depth 
        leggitarget_process(coverage_depth)
}

workflow {
    getSingleInput().cross(getReferences('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }

    module_LEO(input.reads, input.refs)
}

