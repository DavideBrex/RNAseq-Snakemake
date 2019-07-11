rule star:
    input:
        get_trimmed
    output:
        bam   = "02alignments/{sample}/{sample}.bam",
        index = "02alignments/{sample}/{sample}.bam.bai",
        log   = "02alignments/{sample}/Log.final.out"
    log:
       align   = "00log/alignments/{sample}.log",
       rm_dups = "00log/alignments/rm_dup/{sample}.log",
    params:
        out_dir     = "02alignments/{sample}/",
        star_params = config["params"]["star"],
        # path to STAR reference genome index
        index       = config["ref"]["index"]
    threads:
        CLUSTER["star"]["cpu"]
    shadow: 
        "minimal"
    shell: 
        """
        STAR --genomeDir {params.index} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.out_dir} \
        --outSAMtype SAM \
        --outStd SAM \
        {params.star_params} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam} 2>> {log.align}
        """


rule featureCounts:
    input:
         rules.star.output.bam
    output:
        featureCounts = "03featureCounts/{sample}/{sample}.featureCounts",
        summary       = "03featureCounts/{sample}/{sample}.featureCounts.summary",
        counts        = "03featureCounts/{sample}/{sample}.counts",
        rpkm          = "03featureCounts/{sample}/{sample}.rpkm"
        # tpm         = "03featureCounts/{sample}/{sample}.tpm"
    log:
        "00log/featureCounts/{sample}.log"
    params:
        tmp      = "03featureCounts/{sample}/{sample}.seqDepth",
        gtf      = config["ref"]["annotation"],
        options  = config["params"]["featureCounts"],
        # Add -p option for pair-end data if it's the case
        pair_end = lambda w: "-p" if not is_single_end(w.sample) else str()
    threads: 
        CLUSTER["featureCounts"]["cpu"]
    shell:
        """
        featureCounts -T {threads} \
        -a {params.gtf} \
        -o {output.featureCounts} \
        {params.pair_end} {params.options} {input} > {log} 2>&1

        # Create a file with the counts, another for RPKM and another for TPM
        cut -f 1,7 {output.featureCounts} | tail -n +2 > {output.counts}
        seqdepth=$(awk '{{sum+=$2}}END{{print sum}}' {output.counts})
        awk -v s=${{seqdepth}} 'BEGIN{{OFS="\t"}}NR>2{{print $1,(($7)*1000000*1000)/($6*s)}}' {output.featureCounts} > {output.rpkm}
        """
