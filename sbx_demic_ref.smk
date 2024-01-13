try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


DEMIC_FP = MAPPING_FP / "demic"


localrules:
    all_demic_ref,


rule all_demic_ref:
    input:
        expand(DEMIC_FP / "ref_coverage" / "{sample}.txt", sample=Samples.keys()),
        expand(DEMIC_FP / "ref_coverage" / "{sample}.depth", sample=Samples.keys()),
        expand(DEMIC_FP / "ref_coverage" / "{sample}.tsv", sample=Samples.keys()),


rule bowtie2_build_ref:
    input:
        Cfg["sbx_demic_ref"]["refs"],
    output:
        Cfg["sbx_demic_ref"]["refs"] + ".1.bt2",
    conda:
        "envs/sbx_demic_ref_env.yml"
    shell:
        """
        bowtie2-build {input} {input}
        """

rule bowtie2_ref:
    input:
        fasta=Cfg["sbx_demic_ref"]["refs"],
        index=Cfg["sbx_demic_ref"]["refs"] + ".1.bt2",
        reads=expand(
            QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
    output:
        DEMIC_FP / "ref_mapping" / "{sample}.sam",
    conda:
        "envs/sbx_demic_ref_env.yml"
    shell:
        """
        bowtie2 -q -x {input.fasta} -1 {input.reads[0]} -2 {input.reads[1]} -S {output}
        """


rule samtools_sort_ref:
    input:
        DEMIC_FP / "ref_mapping" / "{sample}.sam",
    output:
        sorted_files=DEMIC_FP / "ref_sorted" / "{sample}.sam",
        temp_files=temp(DEMIC_FP / "ref_sorted" / "{sample}-temp.sam")
    threads: 4
    conda:
        "envs/sbx_demic_ref_env.yml"
    shell:
        """
        samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} - -o {output.temp_files}
        samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files}
        """


rule samtools_coverage_ref:
    input:
        DEMIC_FP / "ref_sorted" / "{sample}.sam",
    output:
        hist=DEMIC_FP / "ref_coverage" / "{sample}.txt",
        depth=DEMIC_FP / "ref_coverage" / "{sample}.depth",
        tsv=DEMIC_FP / "ref_coverage" / "{sample}.tsv",
    conda:
        "envs/sbx_demic_ref_env.yml"
    shell:
        """
        samtools coverage {input} -m -o {output.hist}
        samtools coverage {input} -D -o {output.depth}
        samtools coverage {input} -o {output.tsv}
        """