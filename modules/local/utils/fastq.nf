def flowcellLaneFromFastq(path) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab

    // parse first line of a FASTQ file (optionally gzip-compressed)
    // and return the flowcell id and lane number.
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven or eight fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine() // Use the first line in the fastq file for reference
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String machine
    int run_nr
    String fcid
    int lane

    machine = fields[0]
    run_nr = fields[1].toInteger()
    fcid = fields[2]
    lane = fields[3].toInteger()

    [fcid, lane, machine, run_nr]
}


def extractFastqPairFromDir(fastq_path, output){
    // adapted from from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    analysis_id = output.split('/')[-1] // the folder name of params.outdir

    Channel.fromPath("${fastq_path}/*_R1_*.fastq.gz") // Create a channel from all forward reads
        .map{ r1_path ->
            def fastq_files = [r1_path]
            def sample_id = r1_path.getSimpleName().split('_')[0] // extract the sample_id, e.g.,: sample_R1_001.fastq.gz -> sample
            def r2_path = file(r1_path.toString().replace("_R1_", "_R2_")) // reverse reads
            if (r2_path.exists()) { // if reverse reads are available, which should be the case since we sequence paired end
                fastq_files.add(r2_path) // fastq_files -> [r1_path, r2_path]
            } else {
                exit 1, "R2 fastq.gz file not found ${r2_path}"
            }
            def (flowcell, lane) = flowcellLaneFromFastq(r1_path) //extract metadata from reads
            def rg_id = "${sample_id}_${flowcell}_${lane}" // define a readgroup id

            [['id': sample_id, 'rg_id': rg_id, 'flowcell': flowcell, 'analysis_id': analysis_id], fastq_files] // format:  [meta, [r1_path, r2_path]]
        }
}
