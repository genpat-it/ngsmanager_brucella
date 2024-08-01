nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory;extractKey } from '../functions/common'
include { stepInputs;getRisCd } from '../functions/common.nf'
include { getSingleInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def refDir = "/mnt/biowork/databases/REFERENCES/host"

def ex = executionMetadata()

def STEP = '1PP_filtering'
def METHOD = 'bowtie' 

process bowtie2 {
    container "${LOCAL_REGISTRY}/bioinfo/bowtie2:2.1.0--37ad014737"
    containerOptions = "-v ${refDir}:${refDir}:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    when:
      referencePath && referencePath.exists() && !referencePath.empty() && (isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(referencePath)
    output:
      path '*'
      tuple path("${base_ref}.sam"), val(reference), path(referencePath), val(type), emit: sam
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bowtie_${reference}" 
      type = isIlluminaPaired(reads) ? 'paired' : 'single' 
      if (isIlluminaPaired(reads)) {
        """
        bowtie2-build ${referencePath} ${reference}
        bowtie2 -p 8 --very-fast -x ${reference} -1 ${t1} -2 ${t2} -S ${base_ref}.sam 2>> ${base_ref}.log
        """
      } else if (isIonTorrent(reads)) {
        """
        bowtie2-build ${referencePath} ${reference}
        bowtie2 -p 8 --very-fast -x ${reference} -U ${t1} -S ${base_ref}.sam 2>> ${base_ref}.log
        """      
      }     
}

process samtools {
    container 'quay.io/biocontainers/samtools:0.1.19--hf89b575_7'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 18.GB, task.attempt ) }
    input:
      tuple path(sam), val(reference), path(referencePath), val(type)
    output:
      tuple val(riscd), path('*.fastq.gz'), emit: filtered
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base}_R*.fastq.gz"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*_sorted.bam*,*.vcf}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bowtie_${reference}"
      tmpName = "${base}_tmp"
      filteredR1 = "${base}_R1_bowtie_${reference}.fastq.gz"
      filteredR2 = "${base}_R2_bowtie_${reference}.fastq.gz"
      riscd = getRisCd(md, ex, STEP, METHOD)
      if (type == 'paired') {
        """
			  samtools view -bS ${sam} -@ 16 > ${tmpName}.bam 2>> ${base_ref}_samtools.log
		    samtools view -b -f 2 ${tmpName}.bam > ${tmpName}_map.bam 2>> ${base_ref}_samtools.log
		    samtools sort -n ${tmpName}_map.bam ${tmpName}.sorted -@ 16 2>> ${base_ref}_samtools.log
		    samtools bam2fq ${tmpName}.sorted.bam > ${tmpName}.fastq 2>> ${base_ref}_samtools.log
			  grep '^@.*/1\$' -A 3 ${tmpName}.fastq | grep -v -- "^--\$" | gzip > ${filteredR1}
		    grep '^@.*/2\$' -A 3 ${tmpName}.fastq | grep -v -- "^--\$" | gzip > ${filteredR2}   
			  rm ${tmpName}*.*
        """
}     else if (type == 'single') {
        """
		    samtools view -bS ${sam} -@ 16 > ${tmpName}.bam 2>> ${base_ref}_samtools.log
		    samtools view -b -F 4 ${tmpName}.bam > ${tmpName}_map.bam 2>> ${base_ref}_samtools.log
		    samtools sort -n ${tmpName}_map.bam ${tmpName}.sorted -@ 16 2>> ${base_ref}_samtools.log
		    samtools bam2fq ${tmpName}.sorted.bam > ${tmpName}.fastq 2>> ${base_ref}_samtools.log
        gzip -c ${tmpName}.fastq > ${filteredR1} 
        rm ${tmpName}*.*
        """
      }      
}

workflow step_1PP_filtering__bowtie {
    take: 
      reads
      reference 
    main:
      bowtie2(reads, reference) //[ refId, refPath ]
      samtools(bowtie2.out.sam)
    emit:
      samtools.out.filtered      //controlla alla fine
}


workflow {
    getSingleInput().cross(getReference('fa')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }  
    step_1PP_filtering__bowtie(input.reads, input.refs)
}

