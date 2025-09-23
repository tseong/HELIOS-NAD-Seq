#rules/preprocessing.smk

rule run_cutadapt:
    input:
        script="scripts/cutadapt.sh"
    output:
        "cutadapt.done"
    conda: "helios.yml"
    log:   "logs/cutadapt.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_5prime_trim:
    input:
        script="scripts/5prime_trim.sh",
        prev="cutadapt.done"
    output:
        "5prime_trim.done"
    conda: "helios.yml"
    log:   "logs/5prime_trim.log"
    threads: 1
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_trimmomatic:
    input:
        script="scripts/trimmomatic.sh",
        prev="5prime_trim.done"
    output:
        "trimmomatic.done"
    conda: "helios.yml"
    log:   "logs/trimmomatic.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""
