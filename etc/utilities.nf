nextflow.enable.dsl=2

include { filterSampleSheet;sortSampleSheet } from '../functions/samplesheet'


workflow sort_samplesheet {
    sortSampleSheet(params.input, params.output, params.batch_size)
}

workflow filter_samplesheet {
    filterSampleSheet()
}