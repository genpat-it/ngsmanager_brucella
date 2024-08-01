nextflow.enable.dsl=2

include { step_3TX_class__kraken } from '../steps/step_3TX_class__kraken'
include { step_0SQ_rawreads__tenderECDC } from '../steps/step_0SQ_rawreads__tenderECDC'
include { extractKey } from '../functions/common.nf'
include { getDsMetadata } from '../functions/samplesheet.nf'
include { getInput;getSpecies;hasFastqData } from '../functions/parameters.nf'

workflow module_tender_ecdc {
    take: 
        raw_reads
        species
    main:
        raw_reads.branch {
            with_data: hasFastqData(it[1])
            no_reads: true
        }
        .set { rawreads_branched }
        genus_report = step_3TX_class__kraken(rawreads_branched.with_data).genus_report
        rawreads_branched.with_data
            .cross(genus_report) { extractKey(it) }
            .cross(species) { extractKey(it) }
            .multiMap { 
                reads: it[0][0]
                genus_report: it[0][1]
                species: it[1]
            }.set { tender_input }
        step_0SQ_rawreads__tenderECDC(tender_input.reads, tender_input.genus_report, tender_input.species)
    emit:
        no_reads = rawreads_branched.no_reads
        with_data = rawreads_branched.with_data
}

workflow {
    module_tender_ecdc(getInput(), getSpecies())
}    