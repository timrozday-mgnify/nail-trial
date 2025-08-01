/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { CONCATENATE } from  '../modules/local/concatenate/main'
include { FETCHDB } from '../subworkflows/local/fetchdb/main'
include { SEQKIT_TRANSLATE } from '../modules/nf-core/seqkit/translate/main'
include { NAIL_SEARCH } from '../modules/nf-core/nail/search/main'
include { PARSEHMMSEARCHCOVERAGE } from '../modules/local/parsehmmsearchcoverage/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NAILTRIAL {
    main:
    ch_versions = Channel.empty()

    // Fetch databases
    db_ch = Channel
        .from(
            params.databases.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('chunked') && v['chunked']) {
                        v.collect { k_, v_ ->
                            if (v_ instanceof Map) {
                                if (v_.containsKey('base_dir')) {
                                    return [id: k, chunk_id: k_] + v_
                                }
                            }
                        }
                    } else if (v.containsKey('base_dir')) {
                        return [id: k] + v
                    }
                }
            }.flatten()
        )
        .filter { it }

    FETCHDB(db_ch, "${projectDir}/${params.databases.cache_path}")
    dbs_path_ch = FETCHDB.out.dbs

    dbs_path_ch
        .branch { meta, _fp ->
            pfam: meta.id == 'pfam'
        }
        .set { dbs }


    // Parse samplesheet and fetch reads
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

    fasta_ch = samplesheet.map {
        sample, fasta ->
        [
            ['id': sample],
            fasta,
        ]
    }

    SEQKIT_TRANSLATE(fasta_ch)

    pfam_db = dbs.pfam
        .map { meta, fp ->
            file("${fp}/${meta.base_dir}/${meta.files.hmm}")
        }

    ch_chunked_fasta = SEQKIT_TRANSLATE.out.fastx
        .splitFasta(
            size: params.nail_chunksize,
            elem: 1,
            file: true
        )

    ch_chunked_pfam_in = ch_chunked_fasta
        .combine(pfam_db)
        .map { meta, reads, db -> tuple(meta, tuple(reads, db)) }
    
    ch_chunked_pfam_in = ch_chunked_pfam_in
        .groupTuple()
        .flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            def chunksize = chunks_.size()
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunksize), chunk)
            }
        }
        .map { meta, v ->
            def (reads, db) = v
            return [meta, reads, db] 
        }

    NAIL_SEARCH(ch_chunked_pfam_in, false)
    ch_versions = ch_versions.mix(NAIL_SEARCH.out.versions)

    CONCATENATE(
        NAIL_SEARCH.out.target_summary
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}.tbl", results) }
    )
 
    PARSEHMMSEARCHCOVERAGE(
        CONCATENATE.out.concatenated_file
            .combine(pfam_db)
    )

    emit:
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
