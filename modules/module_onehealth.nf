nextflow.enable.dsl=2

include { step_1PP_trimming__fastp } from '../steps/step_1PP_trimming__fastp'
include { step_2AS_denovo__shovill } from '../steps/step_2AS_denovo__shovill'
include { step_3QC_coverage__truecoverage } from '../steps/step_3QC_coverage__truecoverage'
include { step_3TX_species__mash } from '../steps/step_3TX_species__mash'
include { step_3TX_class__confindr } from '../steps/step_3TX_class__confindr'
include { step_4TY_MLST__mlst } from '../steps/step_4TY_MLST__mlst'
include { step_4AN_AMR__resfinder } from '../steps/step_4AN_AMR__resfinder'
include { step_4AN_AMR__abricateefsa } from '../steps/step_4AN_AMR__abricateefsa'
include { step_4TY_cgMLST__chewbbacaefsa } from '../steps/step_4TY_cgMLST__chewbbacaefsa'

include { extractKey;logHeader;isSpeciesSupported } from '../functions/common.nf'
include { getGenusSpecies;getInput;hasFastqData;hasEnoughFastqData } from '../functions/parameters.nf'

def GENUS_SPECIES_ALLOWED = [
  // "Escherichia_coli",
  // "Salmonella_enterica",
  "Listeria_monocytogenes"
]

workflow  module_onehealth {
    take: 
      reads
      genus_species
    main:

      reads.cross(genus_species) { extractKey(it) }
      .filter { isSpeciesSupported(it[1][1], GENUS_SPECIES_ALLOWED, it[0][1], 'module_onehealth') }
      .map { it[0] }
      .set { allowed_reads }

      reads_with_data = allowed_reads.filter { hasFastqData(it[1]) }
      trimmed = step_1PP_trimming__fastp(reads_with_data).trimmed.filter { hasEnoughFastqData(it[1]) }

      step_3TX_species__mash(trimmed)
      assembly = step_2AS_denovo__shovill(trimmed).assembly
      step_4AN_AMR__abricateefsa(assembly)

      // cross each samples with its own genus_species
      trimmed.cross(genus_species) { extractKey(it) }.multiMap { 
        trimmed: it[0]
        genus_species: it[1][1]
      }.set { cross_trimmed }

      step_3QC_coverage__truecoverage(cross_trimmed.trimmed, cross_trimmed.genus_species)
      step_3TX_class__confindr(cross_trimmed.trimmed, cross_trimmed.genus_species)
      step_4AN_AMR__resfinder(cross_trimmed.trimmed, cross_trimmed.genus_species)

      // cross each samples with its own genus_species
      assembly.cross(genus_species) { extractKey(it) }.multiMap { 
        assembly: it[0]
        genus_species: it[1][1]
      }.set { cross_assembly }

      step_4TY_MLST__mlst(cross_assembly.assembly)
      step_4TY_cgMLST__chewbbacaefsa(cross_assembly.assembly, cross_assembly.genus_species, '')
}

workflow {
    module_onehealth(getInput(), getGenusSpecies())
}