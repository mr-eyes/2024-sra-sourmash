import pandas as pd

# Load the list of run IDs from the CSV file
run_ids = pd.read_csv("RunSelector_RunInfo.csv", sep=',', usecols=["Run"])["Run"].tolist()

rule all:
    input:
        expand("SIGS_1K/{sample}.sig", sample=run_ids),
        expand("SIGS_10K/{sample}.sig", sample=run_ids),
        expand("CLEANUP_DONE/cleanup_complete_{sample}.txt", sample=run_ids)

rule sra:
    output:
         "sra/{sample}/{sample}"
    threads: 1
    priority: 0
    resources:
        mem_mb=lambda wildcards, attempt: 5 * 1024 * attempt,
        time=lambda wildcards, attempt: 48 * 60 * attempt,
        runtime=lambda wildcards, attempt: 48 * 60 * attempt,
        partition="bmm"
    shell:
        "aws s3 sync --no-sign-request s3://sra-pub-run-odp/sra/{wildcards.sample} sra/{wildcards.sample}/"

rule fasta:
    input:
        "sra/{sample}/{sample}"
    output:
        "fasta/{sample}.fasta"
    params:
        outdir="fasta/"
    threads: 4
    priority: 10
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        time=lambda wildcards, attempt: 48 * 60 * attempt,
        runtime=lambda wildcards, attempt: 48 * 60 * attempt,
        partition="bmm"
    shell:
        """
        fasterq-dump --skip-technical --fasta-unsorted --threads {threads} --bufsize 1000MB --curcache 10000MB --mem {resources.mem_mb} \
        --outdir {params.outdir} --outfile {wildcards.sample}.fasta.tmp {input}
        mv {params.outdir}{wildcards.sample}.fasta.tmp {params.outdir}{wildcards.sample}.fasta
        """

rule skt1k:
    input:
        "fasta/{sample}.fasta"
    output:
        "SIGS_1K/{sample}.sig"
    params:
        sample_name="{sample}"
    priority: 20
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        time=lambda wildcards, attempt: 16 * 60 * attempt,
        runtime=lambda wildcards, attempt: 16 * 60 * attempt,
        partition="bmm"
    shell:
        "sourmash sketch dna {input} -p k=51,scaled=1000,abund -o {output} --name {params.sample_name}"

rule sk10k:
    input:
        sig="SIGS_1K/{sample}.sig"
    output:
        "SIGS_10K/{sample}.sig"
    params:
        sample_name="{sample}"
    priority: 40
    resources:
        mem_mb=lambda wildcards, attempt: 5 * 1024 * attempt,
        time=lambda wildcards, attempt: 16 * 60 * attempt,
        runtime=lambda wildcards, attempt: 16 * 60 * attempt,
        partition="bmm"
    shell:
        "sourmash signature downsample -q {input.sig} --scaled 10000 -k 51 -o {output}"

rule clean:
    priority: 100
    input:
        sig1="SIGS_1K/{sample}.sig",
        sig2="SIGS_10K/{sample}.sig"
    output:
        "CLEANUP_DONE/cleanup_complete_{sample}.txt"
    shell:
        """
        rm -rf sra/{wildcards.sample}/
        rm -rf fasta/{wildcards.sample}.fasta
        touch {output}
        """