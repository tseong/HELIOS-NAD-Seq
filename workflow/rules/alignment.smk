#rules/alignment.smk

rule run_bowtie2:
    input:
        script="scripts/210425_bowtie2.sh",
        prev="trimmomatic.done"
    output:
        "bowtie2.done"
    conda: "eColiHelios_2.yml"
    log:   "logs/280425_bowtie2.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_featurecounts:
    input:
        script="scripts/280425_featurecounts.sh",
        prev="bowtie2.done"
    output:
        "featurecounts.done"
    conda: "eColiHelios_2.yml"
    log:   "logs/020525_featurecounts.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""
