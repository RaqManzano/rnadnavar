nextflow.enable.dsl = 2

// Note: At the moment it just calls the output file sample.bam and doesn't actually read the sample name
// This assumes that ONT data input is fastq (traditionally ONT outputted fastq, though new basecalling often uses unaligned bam)
// To start from unaligned bam do the following (can change number of threads) - the -T MM,ML allows you to preserve modifications (i.e. 5mC) - this might not be of interest so can skip this but putting here for reference
// samtools fastq -T MM,ML in.fastq | minimap2 -ax map-ont -t 4 ref.fasta - | samtools sort -o out.bam - 

// ONT data might need different trimming and quality filtering to Illumina data

// Variant calling should be ok with GATK (at least for the DNA)

process ONT_DNA_align {
    label 'process_medium'

    conda "bioconda::minimap2"

    input:
    path(DNAfastq)
    path(reference)

    output:
    path "sample_DNA.bam", emit: bam
    path  "ONT_DNA_versions.yml", emit: versions

    script:
    """
    minimap2 -ax map-ont \\
        -t $task.cpus \\
        $reference \\
	    $DNAfastq \\
	| samtools sort -o sample_DNA.bam - 

    samtools index sample_DNA.bam

    cat <<-END_VERSIONS > ONT_DNA_versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
	samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools//g" ) 
   END_VERSIONS
    """

    stub:

    """
    cat <<-END_VERSIONS > ONT_DNA_versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
        samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools//g" )
    END_VERSIONS
    """
}

process ONT_RNA_align {
    label 'process_medium'

    conda "bioconda::minimap2"

    input:
    path(RNAfastq)
    path(reference)

    output:
    path "sample_RNA.bam", emit: bam
    path  "ONT_RNA_versions.yml", emit: versions

    script:
    """
    minimap2 -ax splice \\
        -t $task.cpus \\
        $reference \\
	    $RNAfastq \\
	| samtools sort -o sample_RNA.bam - 

    samtools index sample_RNA.bam

    cat <<-END_VERSIONS > ONT_RNA_versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
	samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools//g" ) 
   END_VERSIONS
    """

    stub:

    """
    cat <<-END_VERSIONS > ONT_RNA_versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
        samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools//g" )
    END_VERSIONS
    """
}

// Paths assume running from the rnadnavar/modules/ONT folder
workflow {

    // To gather used softwares versions 
    versions = Channel.empty()

    // Read in data - test-datasets is in the main rnadnavar repo
    // Pass means these reads passed basecalling with guppy v6.4.6 and sup model (super high accuracy).
    // This data was basecalled with R9.4.1 flow cells which are no longer sold so in the future would be good to update the test dataset to something sequenced with newer ONT chemistry
    // However, hard to find matched samples - these two are matched as far as they are both MCF-7 breast cancer cell line but are not done at the same time or from the same lab etc. 
    DNA_fastq = Channel.fromPath("../../test-datasets/ONT/DNA/pass_chr20.fastq")
    // This test RNA dataset is cDNA - in the future might be good to include direct RNA dataset too but harder to find a public dataset
    RNA_fastq = Channel.fromPath("../../test-datasets/ONT/RNA/SGNex_MCF7_cDNA_replicate1_run3_chr20.fastq")

    // Align RNA to the genome (not transcriptome) hence same reference can be used
    reference = Channel.fromPath("../../test-datasets/ONT/ref/chr20.fa")
    ONT_DNA_align(
        DNA_fastq,
        reference
    )
    ONT_RNA_align(
        RNA_fastq,
        reference
    )

}
