nextflow.enable.dsl=2

// include processes and workflows from other nextflow scripts
include { getInput; getGenusSpecies } from '../functions/parameters.nf'
include { extractDsRef; extractKey } from '../functions/common.nf'
include { module_reads_processing } from '../modules/module_reads_processing.nf'
include { snippy; samtools_pileup; samtools_depth; coverage_minmax } from '../steps/step_2AS_mapping__snippy.nf'
include { confindr } from '../steps/step_3TX_class__confindr.nf'
include { mash_sketch } from '../steps/step_4TY_distance__mash-sketch.nf'


// variables to get reference genome
def referenceCode = 'GCA_000740415.1'
def referencePath = "${params.assets_dir}/module_cansnp_processing/GCF_000740415.1.fasta"


// workflow definition
workflow module_cansnp_processing {
    take: 
      rawReads
      genus_species
    main:
// reads_processing execution (trimming with fastp + QC with fastQC for %Q30)
      module_reads_processing(rawReads)
      trimmedReads = module_reads_processing.out.trimmed_with_data
// snippy execution
      trimmedReads.multiMap {
        trimmed: it
        reference: [ '-', referenceCode, file(referencePath) ]
      }.set { trAndRef }
      bam = snippy(trAndRef.trimmed, trAndRef.reference).bam
// samtools and python scripts execution for quality metrics (Horizontal and Vertical Cov) 
      pileup = samtools_pileup(bam).pileup
      coverage_minmax(bam, 'snippy')
      coverage = samtools_depth(bam, 'snippy').coverage
// confindr execution
      confindr(trimmedReads)
// creation of single-sample msh files using mash sketch
      mash_sketch(trimmedReads)
}


// default workflow (calls created module)
workflow {
    module_cansnp_processing(getInput(), getGenusSpecies())
}