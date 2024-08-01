nextflow.enable.dsl=2

include { step_1PP_trimming__fastp } from '../steps/step_1PP_trimming__fastp'
include { step_2AS_denovo__shovill } from '../steps/step_2AS_denovo__shovill'

include { getSingleInput } from '../functions/parameters.nf'

workflow module_onehealth_assembly {
    take: 
      reads
    main:
      step_1PP_trimming__fastp(reads).trimmed | step_2AS_denovo__shovill
}

workflow {
    module_onehealth_assembly(getSingleInput())
}