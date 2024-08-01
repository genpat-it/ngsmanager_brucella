nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory;stepInputs } from '../functions/common.nf'
include { getSingleInput;getReference;getReferenceCodes;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'snippy' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getExDt(reference, ex) {
    try {        
      reflist = getReferenceCodes()
      if (reflist.size() < 2) {
        return ex.dt
      } 
      if (!reference) {
        log.warn "${reference} should not be empty"
        return ex.dt
      }      
      if (!reflist.contains(reference)) {
        log.warn "${reference} not found in ${reflist}"
        return ex.dt
      }
      return ex.dt + reflist.indexOf(reference)
    } catch(Throwable t) {
        log.warn "getExDt: unexpected exception: ${t.asString()}"
        return ex.dt
    }
}

process snippy {
    container "${LOCAL_REGISTRY}/bioinfo/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    maxForks 4
    when:
      reference_path && reference_path.exists() && !reference_path.empty() && (isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path "**/${base_ref}*"
      tuple path("snippy/${base_ref}.bam"), val(reference), emit: bam
      path '*_input.json'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference, ref_format:ref_format])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{**/*.tab,**/*.fa,**/*.bam*,**/*.bed,**/*.csv,**/*.vcf*,**/*.gff,**/*.html,**/*.txt}', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "*.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      ref_format = (reference_path.getName() ==~ /.+\.f(n)?(a)?(sta)?$/) ? 'fasta' : 'gb'
      is_fasta = r1.getName() ==~ /.+\.fa(sta)?$/
      if (is_fasta) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --ctgs ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """     
      } else if (isIlluminaPaired(reads)) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --R1 ${r1} --R2 ${r2} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """
      } else if (isIonTorrent(reads)) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --se ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """      
      }
}

process samtools_pileup {
    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(bam), val(reference)
    output:
      tuple path("${base_ref}.pileup"), val(reference), emit: pileup    
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_pileup.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      """
      #input should be DS10561267-DT200702_2020.TE.89540.1.2_snippy_NC045512.bam
      samtools mpileup -d 1000 -A -Q 0 ${bam} > ${base_ref}.pileup
      """
}

process samtools_depth {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      tuple path("${base_ref}.coverage"), val(reference), emit: coverage
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_depth.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.coverage', saveAs: { "${base_ref}.coverage.txt" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools depth -a ${bam} | awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' > ${base_ref}.coverage
	    """
}

process coverage_minmax {
    container "${LOCAL_REGISTRY}/bioinfo/samtools:0.1.19--f3869562fe"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      path '*.csv'
      path '*.sh', hidden: true
      tuple path("${base_ref}_samtools_depth.txt"), val(reference), emit: coverage_depth  
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage_minmax.csv"), emit: coverage_extra
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '{*.csv,*.txt}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_minmax.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools view -F 4 -c ${bam} > samtools_view.txt
      samtools depth ${bam} > ${base_ref}_samtools_depth.txt
      /scripts/coverage_minmax.py ${md.cmp} ${md.ds} samtools_view.txt ${base_ref}_samtools_depth.txt ${base_ref}_import_coverage_minmax.csv
	    """
}

process coverage_plot {
    container "quay.io/biocontainers/matplotlib:3.1.2"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      tuple path(coverage_depth), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/result", pattern: '*.png'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_plot.cfg" }
    script:
      md = parseMetadataFromFileName(coverage_depth.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_ivar_${reference}"
      """
      /scripts/coverage_plot.py ${coverage_depth} ${base_ref}_coverage_plot.png
	    """
}

workflow step_2AS_mapping__snippy {
  take: 
    reads
    reference 
  main:
    //snippy(reads, reference)
    bam = snippy(reads, reference).bam
    pileup = samtools_pileup(bam).pileup
    coverage_minmax(bam, 'snippy')
    coverage_minmax.out.coverage_depth | coverage_plot
    coverage = samtools_depth(bam, 'snippy').coverage
  emit:
    //consensus = consensus.map { it[0,2] }
    coverage_depth = coverage_minmax.out.coverage_depth
}

workflow {
    getSingleInput().cross(getReference('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_2AS_mapping__snippy(input.reads, input.refs)
}