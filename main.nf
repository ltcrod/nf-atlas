nextflow.enable.dsl=2
include { SAMTOOLS_VIEW } from './modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SPLITMERGE } from './modules/nf-core/modules/samtools/index/main'
include { ATLAS_RECAL } from './modules/nf-core/modules/atlas/recal/main'
include { ATLAS_SPLITMERGE } from './modules/nf-core/modules/atlas/splitmerge/main'
include { ATLAS_PMD } from './modules/nf-core/modules/atlas/pmd/main'
include { ATLAS_CALL } from './modules/nf-core/modules/atlas/call/main'
include { SAMTOOLS_FAIDX } from './modules/nf-core/modules/samtools/faidx/main'
include { MAKE_PMD_POOLS } from './modules/local/make_pmd_pools'

workflow atlas {

    ch_input_bam = Channel.fromPath( params.input ).map {
        bam -> 
            def meta = [:] //dictionary
            meta.id = bam.baseName.split("\\.").first()
            [meta, bam]
    }
    
    ch_input_read_group = Channel.fromPath( params.read_group ).map {
        txt -> 
            def meta = [:]
            meta.id = txt.baseName.split("\\.").first()
            [meta, txt]
    }

    SAMTOOLS_INDEX ( ch_input_bam )
    
    ch_indexedbam = ch_input_bam
        .join( SAMTOOLS_INDEX.out.bai )
        .join( ch_input_read_group )
        .map {
            meta, bam, bai, read_group ->
                [meta, bam, bai, read_group, []]
        }



    ATLAS_SPLITMERGE ( ch_indexedbam )



    SAMTOOLS_INDEX_SPLITMERGE ( ATLAS_SPLITMERGE.out.data.map {[it[0], it[1]]} )

    //TODO fix issue with piping
    MAKE_PMD_POOLS ( ATLAS_SPLITMERGE.out.data.map {[it[0], it[1]]} )

    //TODO combine splitmerge bam index and RGs

    //ch_bam_indexed2 = ch_input_bam2.merge( SAMTOOLS_INDEX.out.bai )

    //ATLAS_PMD ( ATLAS_SPLITMERGE.out.bam )

}