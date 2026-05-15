include { BWAMEM2_MEM      } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMBAMBA_MARKDUP } from '../../../modules/nf-core/sambamba/markdup/main'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index/main'
include { SEQKIT_SPLIT2    } from '../../../modules/nf-core/seqkit/split2/main'

workflow MAPPING {

    take:
    ch_bwa_index
    ch_genome
    ch_fastq

    main:

    // Split fastq in managable chunks, transform SEQKIT_SPLIT2 output to input per sample fastq part
    SEQKIT_SPLIT2(ch_fastq)
    SEQKIT_SPLIT2.out.reads


    ch_versions = Channel.empty()

    ch_split_fastq = SEQKIT_SPLIT2.out.reads
        .map{ meta, reads -> read_files = reads.sort(false) {
                a,b -> a.getName()
                    .tokenize('.')[1] <=> b.getName()
                    .tokenize('.')[1]
                }
             .collate(2)
             [meta, read_files]
        }
        .transpose()
        .map{ meta, read_files -> [meta + [
                    split_fastq_part:read_files[0]
                        .getName()
                        .tokenize('.')[1]
                ], read_files
            ]
        }

    // Mapping
    BWAMEM2_MEM(
        ch_split_fastq,
        ch_bwa_index,
        ch_genome,
        true
    )

    SAMBAMBA_MARKDUP(
        BWAMEM2_MEM.out.bam
            .map{ meta, bam -> [meta - meta.subMap('rg_id', 'flowcell', 'split_fastq_part'), bam]}
            .groupTuple()
    )
    SAMTOOLS_INDEX(
        SAMBAMBA_MARKDUP.out.bam
    )

    ch_versions = channel.versions.mix(SEQKIT_SPLIT2.out.versions)
    ch_versions = channel.versions.mix(BWAMEM2_MEM.out.versions)
    ch_versions = channel.versions.mix(SAMBAMBA_MARKDUP.out.versions)
    ch_versions = channel.versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam = SAMBAMBA_MARKDUP.out.bam
    bai = SAMTOOLS_INDEX.out.bai
    sambamba_txt = SAMBAMBA_MARKDUP.out.txt
    versions = ch_versions
}
