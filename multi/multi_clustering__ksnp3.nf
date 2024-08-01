nextflow.enable.dsl=2

include { getInput;param } from '../functions/parameters.nf'
include { taskMemory;flattenPath } from '../functions/common.nf'

def GEO_RESOLUTION_COLUMNS='Comune'
def DATE_COLUMN_1='Sampling Date'
def DATE_COLUMN_2='DataPrelievo'
def SAMPLE_COLUMN='CMP'

process ksnp3 {
    container "${LOCAL_REGISTRY}/bioinfo/ksnp3:3.0--addd2c2d0e"
    input:
      path(assembly)
      val(kmers_size)
      val(analysis_type)
    output:
      path("results/${outfile}"), emit: fasta
      path '*.sh', hidden: true
    publishDir mode: 'copy', "${params.outdir}", pattern: "*results/${outfile}", saveAs: { "matrix.fasta" }   
    publishDir mode: 'copy', "${params.outdir}", pattern: '*results/*.fasta', saveAs: { filename -> flattenPath(filename) }   
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "ksnp3.cfg" }
    script:
      extra_params = (analysis_type == 'core' ? '-core' : '')
      outfile = (analysis_type == 'core' ? 'core_SNPs_matrix.fasta' : 'SNPs_all_matrix.fasta')
      """
        trap "rm -Rf results/TemporaryFilesToDelete" EXIT
        for FILE in ${assembly}; do \
          CMP=`echo \$FILE | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/' | sed 's/\\./-/g'`; \
          ln -s \$FILE \${CMP}.fasta ; \
          echo -e "`pwd`/\${CMP}.fasta\t\$CMP" >> input.tsv ; \
        done        
        kSNP3 -in input.tsv -k ${kmers_size} -NJ ${extra_params} -outdir results > ksnp3.log
      """
}

process iqtree {
    container "quay.io/biocontainers/iqtree:1.6.12--he513fc3_0"
    input:
      path(fasta)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path("*.treefile"), emit: nwk
    publishDir mode: 'copy', "${params.outdir}", pattern: '*.treefile', saveAs: { "matrix.nwk" } 
    publishDir mode: 'copy', "${params.outdir}", pattern: '*.treefile', saveAs: { filename -> filename.replace(".treefile", ".nwk") } 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "iqtree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "iqtree.cfg" }
    script:
      """
        iqtree -nt AUTO -s ${fasta}      
        sed -i 's/-/./g' *.treefile
      """
}

process grapetree {
    container "quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(matrix)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '*.nwk'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "grapetree.cfg" }
    script:
      """
        grapetree -p ${matrix} > matrix.nwk
      """
}

process augur {
    container "quay.io/biocontainers/augur:22.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path(nwk)
      path(metadata)
      path(geodata)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur.cfg" }
    script:
      """
        cat ${metadata} | sed 's/${SAMPLE_COLUMN}/name/i' \
           | sed 's/${DATE_COLUMN_1}/date/i' \
           | sed 's/${DATE_COLUMN_2}/date/i' > augur_metadata.tsv
        METADATA_LIST=\$(head -n 1 augur_metadata.tsv | tr \$'\t' ' ')
        augur refine --tree ${nwk} --output-tree tree_tt.nwk --output-node-data refine.node.json --metadata augur_metadata.tsv
        augur export v2 --tree tree_tt.nwk --node-data refine.node.json --output auspice.json \
          --color-by-metadata \${METADATA_LIST} \
          --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
          --metadata augur_metadata.tsv \
          --lat-longs ${geodata}
      """
}

workflow multi_clustering__ksnp3 {
    take: 
        input
        kmers_size
        analysis_type
        metadata
        geodata
    main:
        input
          .map { it[1] }
          .collect()
          .set { inputSet }
        matrix = ksnp3(inputSet, kmers_size, analysis_type).fasta
        iqtree(matrix)
        grapetree(matrix)
        augur(iqtree.out.nwk, metadata, geodata)
}

workflow {
    multi_clustering__ksnp3(getInput(), param('kmers_size'), param('analysis_type'), param('metadata'), param('geodata'));
}
