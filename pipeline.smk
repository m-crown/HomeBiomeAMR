rule all:
    input:
        expand( config["outdir"] + "/human_removed/{sample}_clean_{read}.fastq.gz", sample=config["samples"], read = ["R1", "R2"]),
        expand( config["outdir"] + "/shortBRED/{sample}/{sample}_shortbred_results.tsv", sample = config["samples"]),
        expand( config["outdir"] + "/kraken/{sample}_{level}.bracken", sample = config["samples"], level=config["levels"]),
        expand( config["outdir"] + "/kraken/merged_bracken_{level}.txt", level = config["levels"]),
        expand(config["outdir"] + "/kraken/bracken_plot_{level}.png", level = config["levels"]),
        expand( config["outdir"] + "/kraken/{sample}.breport", sample = config["samples"]),
        config["outdir"] + "/read_stats/read_stats_plot.png",
        resistance_samples = config["outdir"] + "/shortBRED_analysis/all_samples_resistance_mechanism.html",
         resistance_combined = config["outdir"] + "/shortBRED_analysis/combined_resistance_mechanism.html",
         drug_samples = config["outdir"] + "/shortBRED_analysis/all_samples_drug_resistance.html",
         drug_combined = config["outdir"] + "/shortBRED_analysis/combined_drug_resistance.html"

rule trim_galore:
    input:
        r1= config["indir"] + "/{sample}_ME_L001_R1_001.fastq.gz",
        r2= config["indir"] + "/{sample}_ME_L001_R2_001.fastq.gz"
    output:
        trimmed_r1= config["outdir"] + "/trim_galore/{sample}_val_1.fq.gz",
        trimmed_r2= config["outdir"] + "/trim_galore/{sample}_val_2.fq.gz"
    params:
        length=100,
        trimgalore = config["trimgalore"],
        out = config["outdir"] + "/trim_galore"
    threads: 4
    shell:
        "{params.trimgalore} --basename {wildcards.sample} --cores {threads} --length {params.length} --paired -o {params.out} {input.r1} {input.r2}"

rule bbmap:
    input:
        r1= config["outdir"] + "/trim_galore/{sample}_val_1.fq.gz",
        r2= config["outdir"] + "/trim_galore/{sample}_val_2.fq.gz"
    output:
        clean= config["outdir"] + "/human_removed/{sample}_clean.fq.gz",
        human= config["outdir"] + "/human_removed/{sample}_human.fq.gz"
    params:
        bbmap_dir = config["bbmap_dir"],
        remove_human_dir = config["remove_human"]
    threads: 20
    shell:
        "{params.bbmap_dir}/bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path={params.remove_human_dir} qtrim=rl trimq=10 untrim -Xmx23g t={threads} in1={input.r1} in2={input.r2} outu={output.clean} outm={output.human}"

rule reformat:
    input:
        clean= config["outdir"] + "/human_removed/{sample}_clean.fq.gz"
    output:
        r1= config["outdir"] + "/human_removed/{sample}_clean_R1.fastq.gz",
        r2= config["outdir"] + "/human_removed/{sample}_clean_R2.fastq.gz"
    params:
        bbmap_dir = config["bbmap_dir"]
    shell:
        "{params.bbmap_dir}/reformat.sh in={input.clean} out1={output.r1} out2={output.r2}"
        
rule kraken:
    input:
        r1= config["outdir"] + "/human_removed/{sample}_clean_R1.fastq.gz",
        r2= config["outdir"] + "/human_removed/{sample}_clean_R2.fastq.gz"
    output:
        kraken_file = config["outdir"] + "/kraken/{sample}_classified.kraken",
        kraken_report = config["outdir"] + "/kraken/{sample}_classified.kreport",
        classified_reads_1 = config["outdir"] + "/kraken/{sample}_classified_1.fastq",
        classified_reads_2 = config["outdir"] + "/kraken/{sample}_classified_2.fastq"
    params:
        kraken_db = config["kraken_db"]
    threads: 40
    shell:
        "kraken2 --db {params.kraken_db} --threads {threads} --paired --output {output.kraken_file} --report {output.kraken_report} --classified-out {config[outdir]}/kraken/{wildcards.sample}_classified#.fastq {input.r1} {input.r2}"

rule select_bacterial_reads:
    input:
        kraken_file = config["outdir"] + "/kraken/{sample}_classified.kraken",
        kraken_report = config["outdir"] + "/kraken/{sample}_classified.kreport",
        classified_reads_1 = config["outdir"] + "/kraken/{sample}_classified_1.fastq",
        classified_reads_2 = config["outdir"] + "/kraken/{sample}_classified_2.fastq"
    output:
        bacterial_reads_1 = config["outdir"] + "/kraken/{sample}_bacterial_1.fasta",
        bacterial_reads_2 = config["outdir"] + "/kraken/{sample}_bacterial_2.fasta"
    params:
    shell:
        """python3 extract_kraken_reads.py --include-children -t 2 -k {input.kraken_file} -r {input.kraken_report} -s1 {input.classified_reads_1} -s2 {input.classified_reads_2} -o {output.bacterial_reads_1} -o2 {output.bacterial_reads_2}"""
        
rule kraken_bacterial:
    input:
        r1= config["outdir"] + "/kraken/{sample}_bacterial_1.fasta",
        r2= config["outdir"] + "/kraken/{sample}_bacterial_2.fasta"
    output:
        kraken_file = config["outdir"] + "/kraken/{sample}.kraken",
        kraken_report = config["outdir"] + "/kraken/{sample}.kreport",
    params:
        kraken_db = config["kraken_db"]
    threads: 40
    shell:
        "kraken2 --db {params.kraken_db} --threads {threads} --paired --output {output.kraken_file} --report {output.kraken_report} {input.r1} {input.r2}"

rule bracken:
    input:
        kraken_report = config["outdir"] + "/kraken/{sample}.kreport"
    output:
        bracken_file= config["outdir"] + "/kraken/{sample}_{level}.bracken",
        bracken_report = config["outdir"] + "/kraken/{sample}_{level}.breport"
    params:
        kraken_db = config["kraken_db"],
        read_length = config["read_length"],
        bracken_level = config["bracken_level"],
        levels=config["levels"]
    shell:
        """bracken -d {params.kraken_db} -i {input.kraken_report} -r {params.read_length} -l {params.bracken_level} -l {wildcards.level} -t {threads} -o {output.bracken_file} -w {output.bracken_report}"""

rule combine_bracken_outputs:
    input:
        bracken_files=expand(config["outdir"] + "/kraken/{sample}_{level}.bracken", sample=config["samples"], level=config["levels"])
    output:
        config["outdir"] + "/kraken/merged_bracken_{level}.txt"
    params:
        dir = config["outdir"] + "/kraken",
        names = ','.join(sorted(config["samples"]))
    shell:
        """
        python3 combine_bracken_outputs.py --files {params.dir}/*_{wildcards.level}.bracken --names {params.names} --output {output}
        """
rule plot_bracken_abundances:
    input:
        expand(config["outdir"] + "/kraken/merged_bracken_{level}.txt", level = config["levels"])
    output:
        expand(config["outdir"] + "/kraken/bracken_plot_{level}.png", level = config["levels"])
    params:
        dir = config["outdir"] + "/kraken/"
    shell:
        """Rscript bracken_plot.R --dir {params.dir}"""

rule shortbred_quantify:
    input:
        r1= config["outdir"] + "/kraken/{sample}_bacterial_1.fasta",
        r2 = config["outdir"] + "/kraken/{sample}_bacterial_2.fasta"
    output:
        results = config["outdir"] + "/shortBRED/{sample}/{sample}_shortbred_results.tsv"
    params:
        markers= config["shortbred_markers"],
        tmp= config["outdir"] + "/shortBRED/{sample}/"
    threads: 20
    shell:
        "shortbred_quantify.py --markers {params.markers} --wgs {input.r1} {input.r2} --results {output.results} --tmp {params.tmp} --threads {threads}"

rule shortbred_analyse:
    input:
        expand( config["outdir"] + "/shortBRED/{sample}/{sample}_shortbred_results.tsv", sample=config["samples"])
    output:
        drug_file = expand( config["outdir"] + "/shortBRED_analysis/{sample}_drug.tsv", sample=config["samples"]) ,
        resistance_file = expand(config["outdir"] + "/shortBRED_analysis/{sample}_res.tsv", sample=config["samples"]),
        combined_drug = config["outdir"] + "/shortBRED_analysis/rpkms_sum_drug_combined.tsv",
        combined_resistance = config["outdir"] + "/shortBRED_analysis/rpkms_sum_res_combined.tsv"
    params:
        metadata = config["metadata"],
        shortbred_dir = config["outdir"] + "/shortBRED",
        out_dir = config["outdir"] + "/shortBRED_analysis/",
        aro_index = config["card_aro_index"]
        
    shell:
        "Rscript shortbred_analysis_2023.R --metadata {params.metadata} --shortbred_dir {params.shortbred_dir} --out_dir {params.out_dir} --aro_index {params.aro_index}"

rule krona_plots:
    input:
         samples_drug = expand( config["outdir"] + "/shortBRED_analysis/{sample}_drug.tsv", sample=config["samples"]),
         samples_resistance = expand( config["outdir"] + "/shortBRED_analysis/{sample}_res.tsv", sample=config["samples"]),
         combined_drug = config["outdir"] + "/shortBRED_analysis/rpkms_sum_drug_combined.tsv",
         combined_resistance = config["outdir"] + "/shortBRED_analysis/rpkms_sum_res_combined.tsv"
    output:
         resistance_samples = config["outdir"] + "/shortBRED_analysis/all_samples_resistance_mechanism.html",
         resistance_combined = config["outdir"] + "/shortBRED_analysis/combined_resistance_mechanism.html",
         drug_samples = config["outdir"] + "/shortBRED_analysis/all_samples_drug_resistance.html",
         drug_combined = config["outdir"] + "/shortBRED_analysis/combined_drug_resistance.html"
    shell:
         """ktImportText -o {output.resistance_samples} {input.samples_resistance}; ktImportText -o {output.drug_samples} {input.samples_drug}; ktImportText -o {output.resistance_combined} {input.combined_resistance}; ktImportText -o {output.drug_combined} {input.combined_drug}"""
         
rule read_stats:
    input:
        bacterial_reads = expand(config["outdir"] + "/kraken/{sample}_bacterial_{read_num}.fasta", sample = config["samples"], read_num = [1,2])
    output:
        raw_read_stats = config["outdir"] + "/read_stats/raw_read_stats.tsv",
        trim_galore_read_stats = config["outdir"] + "/read_stats/trim_galore_read_stats.tsv",
        bbmap_read_stats = config["outdir"] + "/read_stats/bbmap_read_stats.tsv",
        kraken_stats = config["outdir"] + "/read_stats/kraken_read_stats.tsv"
    params:
        trim_galore_dir = config["outdir"] + "/trim_galore",
        bbmap_dir = config["outdir"] + "/human_removed",
        raw_dir = config["indir"],
        kraken_dir = config["outdir"] + "/kraken"
    threads: 20
    shell:
        """seqkit stats {params.raw_dir}/*.fastq.gz -T -j {threads} -o {output.raw_read_stats}; seqkit stats {params.trim_galore_dir}/*.fq.gz -T -j {threads} -o {output.trim_galore_read_stats}; seqkit stats {params.bbmap_dir}/*.fastq.gz -T -j {threads} -o {output.bbmap_read_stats}; seqkit stats {input.bacterial_reads} -T -j {threads} -o {output.kraken_stats}"""
        
rule read_plot:
    input:
        raw_read_stats = config["outdir"] + "/read_stats/raw_read_stats.tsv",
        trim_galore_read_stats = config["outdir"] + "/read_stats/trim_galore_read_stats.tsv",
        bbmap_read_stats = config["outdir"] + "/read_stats/bbmap_read_stats.tsv",
        kraken_stats = config["outdir"] + "/read_stats/kraken_read_stats.tsv"
    output:
        config["outdir"] + "/read_stats/read_stats_plot.png"
    shell:
        """Rscript read_distributions.R"""