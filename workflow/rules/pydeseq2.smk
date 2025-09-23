# rules/pydeseq2.smk

TP_LIST = [f"tp{i}" for i in range(1,17)]

rule merge_featurecounts:
    input:
        script="scripts/merge_featureCounts_barcode.py",
        prev="featurecounts.done"
    output:
        temp("merge_featurecounts_{tp}.done"),
        counts="table_Astart/{tp}/merged_by_barcode_Astart_readCount.csv"
    conda: "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/eColiHelios_2.yml"
    log:   "logs/140825_merge_featurecounts_{tp}.log"
    threads: 1
    shell: r"""python {input.script} > {log} 2>&1; touch {output[0]}"""

rule filter_dds:
    input:
        script="scripts/060525_filter_dds.py",
        counts="table_Astart/{tp}/merged_by_barcode_Astart_readCount.csv",
        prev="merge_featurecounts_{tp}.done"
    output:
        temp("filter_dds_{tp}.done")
    conda: "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/eColiHelios_2.yml"
    log:   "logs/140825_filter_dds_{tp}.log"
    threads: 1
    shell: r"""python {input.script} > {log} 2>&1; touch {output}"""

rule stat_test:
    input:
        script="scripts/060525_stat_test.py",
        counts="table_Astart/{tp}/merged_by_barcode_Astart_readCount.csv",
        prev="filter_dds_{tp}.done"
    output:
        "table_Astart/{tp}/pydeseq2_results_merged_by_barcode_Astart_readCount.csv"
    conda: "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/eColiHelios_2.yml"
    log:   "logs/140825_stat_test_{tp}.log"
    threads: 1
    shell: r"""python {input.script} > {log} 2>&1"""
