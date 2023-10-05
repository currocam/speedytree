import random
import pandas as pd
import numpy as np

DATA_DIR = "pfam_alignments"
(SAMPLE,) = glob_wildcards(f"{DATA_DIR}/{{sample}}.sth")


rule all:
    input:
        expand(
            "benchmarks/speedytree/{sample}{seed}_{algo}_{threads}.txt",
            algo=["naive", "rapidnj", "hybrid"],
            sample=SAMPLE,
            seed=["", "_seed_1", "_seed_2"],
            threads=[1, 8],
        ),
        expand(
            "benchmarks/{program}/{sample}{seed}.txt",
            program=["rapidnj", "quicktree"],
            sample=SAMPLE,
            seed=["", "_seed_1", "_seed_2"],
        ),


rule speedytree_binary:
    output:
        "target/release/speedytree",
    shell:
        "cargo build --release"


rule speedytree:
    input:
        "target/snakemake/{sample}.phy",
        "target/release/speedytree",
    benchmark:
        "benchmarks/speedytree/{sample}_{algo}_{threads}.txt"
    shell:
        "./{input[1]} --{wildcards.algo} -c {wildcards.threads} < {input[0]} > /dev/null"


rule rapidnj:
    input:
        "target/snakemake/{sample}.phy",
    benchmark:
        "benchmarks/rapidnj/{sample}.txt"
    shell:
        "rapidnj -i pd {input} > /dev/null"


rule quicktree:
    input:
        "target/snakemake/{sample}.phy",
    benchmark:
        "benchmarks/quicktree/{sample}.txt"
    shell:
        "quicktree -in m {input} > /dev/null"


rule phylip:
    input:
        f"{DATA_DIR}/{{sample}}.sth",
    output:
        "target/snakemake/{sample}.phy",
    shell:
        "rapidnj -i sth -o m {input} > {output}"


rule permutate_phylip:
    input:
        "target/snakemake/{sample}.phy",
    output:
        "target/snakemake/{sample}_seed_{seed}.phy",
    run:
        random.seed(wildcards.seed)
        data = pd.read_table(input[0], delim_whitespace=True, skiprows=1, header=None)
        data.set_index(0, inplace=True)
        names = data.index
        data = data.values
        # Set column names as index
        n = data.shape[1]
        perm = list(range(n))
        random.shuffle(perm)
        perm_data = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                perm_data[i, j] = data[perm[i], perm[j]]
        perm_names = [names[i] for i in perm]
        with open(output[0], "w") as f:
            f.write(f"{perm_data.shape[0]}\n")
            for name, row in zip(perm_names, perm_data):
                name = name.strip()
                f.write(f"{name} {' '.join(map(str, row))}\n")
