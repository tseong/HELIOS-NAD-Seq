#rules/alignment.smk

rule run_bowtie2:
    input:
        script="scripts/bowtie2.sh",
        prev="trimmomatic.done"
    output:
        "bowtie2.done"
    conda: "helios.yml"
    log:   "logs/bowtie2.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_featurecounts:
    input:
        script="scripts/featurecounts.sh",
        prev="bowtie2.done"
    output:
        "featurecounts.done"
    conda: "helios.yml"
    log:   "logs/featurecounts.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""
