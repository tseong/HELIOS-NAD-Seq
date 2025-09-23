#rules/preprocessing.smk

rule run_cutadapt:
    input:
        script="scripts/230425_cutadapt.sh"
    output:
        "cutadapt.done"
    conda: "eColiHelios_2.yml"
    log:   "logs/180425_cutadapt.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_5prime_trim:
    input:
        script="scripts/200425_5prime_trim.sh",
        prev="cutadapt.done"
    output:
        "5prime_trim.done"
    conda: "eColiHelios_2.yml"
    log:   "logs/230425_5prime_trim.log"
    threads: 1
    shell: r"""bash {input.script} > {log} 2>&1"""

rule run_trimmomatic:
    input:
        script="scripts/200425_trimmomatic.sh",
        prev="5prime_trim.done"
    output:
        "trimmomatic.done"
    conda: "eColiHelios_2.yml"
    log:   "logs/240425_trimmomatic.log"
    threads: 4
    shell: r"""bash {input.script} > {log} 2>&1"""
