import numpy as np

configfile: "config/config.yaml"
url_clumpak = config['url_clumpak']
url_sgdp = config['url_sgdp']
url_metadata_sgdp = config['url_metadata_sgdp']
url_hg19_to_hg38_chain = config["url_hg19_to_hg38_chain"]
url_hg38_fasta = config['url_hg38_fasta']
url_hollfelder = config['url_hollfelder']
url_arauna_bim = config['url_arauna_bim']
url_arauna_bed = config['url_arauna_bed']
url_arauna_fam = config['url_arauna_fam']
geo_coords_arauna = config['geo_coords_arauna']
populations_of_interest_sgdp = config['populations_of_interest_sgdp']
populations_to_exclude_hollfelder = config['populations_to_exclude_hollfelder']
data_path =  config["data_path"]
results_path = config['results_path']
plink_prefix_crawford_and_scheinfeldt = config['plink_prefix_crawford_and_scheinfeldt']
geo_coords_crawford_and_scheinfeldt = config['geo_coords_crawford_and_scheinfeldt']
geo_coords_hollfelder = config['geo_coords_hollfelder']
plink_prefix_fortes_lima = config['plink_prefix_fortes_lima']
geo_coords_fortes_lima = config['geo_coords_fortes_lima']
samples_downsampling = config['samples_downsampling']
maf = config['maf'] # minor allele frequency
geno = config['geno'] # max missing genotype rate per variant
max_r2 = config['max_r2']
K = np.arange(config['min_K'], config['max_K'] + 1)
bcftools_path = config['bcftools_path']

chromosomes = np.arange(1, 23)

rule all:
    input:
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_pca.png",
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_feems_plot.pdf",
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_plots.pdf",
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_plot_best_K.pdf",
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_kriging_best_K.pdf"

rule download_and_install_clumpak:
    output:
        directory(url_clumpak.split('/')[-1].split('.zip')[0])
    params:
        url = url_clumpak,
        zip = url_clumpak.split('/')[-1]
    shell:
        "wget -q {params.url}; unzip {params.zip}; rm {params.zip}; mv {output} {output}_tmp; "
        "unzip -d {output}_tmp {output}_tmp/26_03_2015_CLUMPAK.zip; mv {output}_tmp/26_03_2015_CLUMPAK/CLUMPAK .;"
        "rm -r {output}_tmp; chmod +x {output}/distruct/*; chmod +x {output}/CLUMPP/*; chmod +x {output}/mcl/*"

rule download_sgdp_vcfs:
    output:
        temp(data_path + "SGDP/" + url_sgdp.split('/')[-1])
    params:
        url = url_sgdp,
        path = data_path + "SGDP/"
    shell:
        "wget -q -P {params.path} {params.url}"

rule download_sgdp_metadata:
    output:
        data_path + "SGDP/" + url_metadata_sgdp.split('/')[-1]
    params:
        url = url_metadata_sgdp,
        path = data_path + "SGDP/",
        pops_of_interest = "|".join(populations_of_interest_sgdp)
    shell:
        "wget -q -O- {params.url} | awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{print $3, $5, $6, $7, $8, $12, $13, $2}}' | "
        "grep -E \"{params.pops_of_interest}\" > {output}"

checkpoint extract_sgdp_vcf:
    input:
        data_path + "SGDP/" +  url_sgdp.split('/')[-1]
    output:
        directory(data_path + "SGDP/vcfs/")
    params:
        dir = data_path + "SGDP/vcfs/"
    shell:
        "mkdir -p  {params.dir}; tar xvf {input} -C {params.dir}"

def get_input_merge_sgdp_vcfs(wildcards):
    data_dir = checkpoints.extract_sgdp_vcf.get(**wildcards).output[0]
    samples = sorted(glob_wildcards(data_dir + '/{sample}.vcf.gz').sample)
    vcfs = [data_dir + f'/{sample}.vcf.gz' for sample in samples]
    return vcfs

rule merge_sgdp_vcfs:
    input:
        get_input_merge_sgdp_vcfs
    output:
        data_path + 'SGDP/ALL.vcf.gz'
    params:
        bcftools = bcftools_path,
        chrom = ",".join(chromosomes.astype(str).tolist())
    threads: 64
    shell:
        "{params.bcftools} merge -0 -r {params.chrom} -Oz -o {output} --threads {threads} {input}"

rule convert_vcf_to_plink_sgdp:
    input:
        ancient(data_path + 'SGDP/ALL.vcf.gz')
    output:
        multiext(data_path + 'SGDP/ALL', ".bed", ".bim", ".fam")
    conda:
        "envs/plink.yaml"
    params:
        base = data_path + 'SGDP/ALL'
    threads: 64
    shell:
        "plink2 --vcf {input} --snps-only --max-alleles 2 --make-bed --set-all-var-ids @:# "
        "--out {params.base} --threads {threads}"

rule select_samples_sgdp:
    input:
        fam = data_path + "SGDP/ALL.fam",
        metadata =data_path + "SGDP/" + url_metadata_sgdp.split('/')[-1]
    output:
        temp("to_keep_sgdp.tab")
    shell:
        "cut -f2 {input.fam} | tr \"\n\" \"|\" | xargs -I {{}} grep -Ew {{}} {input.metadata} | "
        "awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{print 0, $NF}}' > {output};"


rule extract_samples_sgdp:
    input:
        bed = data_path + "SGDP/ALL.bed",
        samples = 'to_keep_sgdp.tab'
    output:
        temp(multiext(data_path + "SGDP/ALL_samples_of_interest", ".bed", ".bim", ".fam"))
    params:
        input_base = data_path + "SGDP/ALL",
        output_base = data_path + "SGDP/ALL_samples_of_interest"
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        "plink2 --bfile {params.input_base} --keep {input.samples} --make-bed --out {params.output_base} --threads {threads}"

rule get_updated_ids_sgdp:
    input:
        old = 'to_keep_sgdp.tab',
        metadata =data_path + "SGDP/" + url_metadata_sgdp.split('/')[-1]
    output:
        temp("updated_sample_ids_sgdp.tab")
    shell:
        "cut -f2 {input.old} | tr \"\n\" \"|\" | xargs -I {{}} grep -E {{}} {input.metadata} | "
        "awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{print $3, $2}}' | paste {input.old} - > {output}"

rule update_sample_ids_sgdp:
    input:
        ids="updated_sample_ids_sgdp.tab",
        bed=data_path + "SGDP/ALL_samples_of_interest.bed",
        fam=data_path + "SGDP/ALL_samples_of_interest.fam",
        bim=data_path + "SGDP/ALL_samples_of_interest.bim"
    output:
        multiext(data_path + "SGDP/ALL_samples_of_interest_updated_ids", ".bed", ".bim", ".fam")
    params:
        input_base=data_path + "SGDP/ALL_samples_of_interest",
        output_base=data_path + "SGDP/ALL_samples_of_interest_updated_ids"
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        'plink2 --bfile {params.input_base} --update-ids {input.ids} --make-bed --out {params.output_base} '
        '--threads {threads}'

rule get_samples_for_downsampling_sgdp:
    input:
        fam=data_path + "SGDP/ALL_samples_of_interest_updated_ids.fam"
    output:
        temp("downsampled_sgdp.tab")
    params:
        n_samples = samples_downsampling
    shell:
        "cut -f1 {input.fam} | sort | uniq | xargs -I {{}} bash -c \"grep {{}} {input.fam} | cut -f1,2 | "
        "shuf -n{params.n_samples}\" bash > {output}"

rule downsample_sgdp:
    input:
        "downsampled_sgdp.tab",
        multiext(data_path + "SGDP/ALL_samples_of_interest_updated_ids", ".bed", ".fam", ".bim")
    output:
        multiext(data_path + "SGDP/downsampled_updated_ids",".bed",".fam",".bim")
    params:
        input_base = data_path + "SGDP/ALL_samples_of_interest_updated_ids",
        output_base = data_path + "SGDP/downsampled_updated_ids"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --keep {input[0]} --make-bed --out {params.output_base} --threads {threads}"


rule get_samples_with_geo_data_available_crawford_and_scheinfeldt:
    input:
        coords=data_path + geo_coords_crawford_and_scheinfeldt
    output:
        to_keep=temp("crawford_and_scheinfeldt_geo_data_avail.tab")
    shell:
        "cut -f1,2 {input.coords} > {output.to_keep}"

rule remove_samples_without_geo_data_crawford_and_scheinfeldt:
    input:
        keep="crawford_and_scheinfeldt_geo_data_avail.tab",
        fam=data_path + plink_prefix_crawford_and_scheinfeldt + ".fam"
    output:
        temp(multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_geo_avail",".bed",".fam",".bim"))
    params:
        input_base = data_path + plink_prefix_crawford_and_scheinfeldt,
        output_base = data_path + plink_prefix_crawford_and_scheinfeldt + "_geo_avail"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --keep {input.keep} --make-bed --out {params.output_base} --threads {threads}"

rule get_updated_fids_crawford_and_scheinfeldt:
    input:
        fam=data_path + plink_prefix_crawford_and_scheinfeldt + "_geo_avail.fam",
        meta=data_path + geo_coords_crawford_and_scheinfeldt
    output:
        temp("updated_fids_crawford_and_scheinfeldt.tab")
    shell:
        "cut -f2 {input.fam} | xargs -I {{}} grep {{}} {input.meta} | "
        "awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{print $3, $2}}' | sed -e 's/\///g' | "
        "paste <(cut -f1,2 {input.fam} | sed 's/\\s/\\t/g') - > {output}"

rule update_fids_crawford_and_scheinfeldt:
    input:
        ids="updated_fids_crawford_and_scheinfeldt.tab",
        fam=multiext(data_path+ plink_prefix_crawford_and_scheinfeldt + "_geo_avail", ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_updated_fids", ".bed", ".fam", ".bim"))
    params:
        input_base=data_path + plink_prefix_crawford_and_scheinfeldt + "_geo_avail",
        output_base=data_path + plink_prefix_crawford_and_scheinfeldt + "_updated_fids"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --update-ids {input.ids} --make-bed --threads {threads} "
        "--out {params.output_base}"

rule select_samples_crawford_and_scheinfeldt:
    input:
        samples = "updated_fids_crawford_and_scheinfeldt.tab"
    output:
        temp(data_path + "crawford_and_scheinfeldt_sampled.tab")
    params:
        samples = "updated_fids_crawford_and_scheinfeldt.tab",
        n_samples=samples_downsampling
    shell:
        "cat {input.samples} | sed -e 's/\///g' -e 's/|/|\\\/' | "
        "cut -f3 | sort | uniq | xargs -I {{}} bash -c \"grep -w \\\"{{}}\\\" <( cat {input.samples} | sed -e 's/\///g') | "
        "shuf -n {params.n_samples} | cut -f3,4 >> {output}\" sh"

rule downsample_crawford_and_scheinfeldt:
    input:
        samples=data_path + "crawford_and_scheinfeldt_sampled.tab",
        bed=multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_updated_fids", ".bed", ".bim", ".fam")
    output:
        multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled", ".bed", ".bim", ".fam")
    params:
        input_base = data_path + plink_prefix_crawford_and_scheinfeldt + "_updated_fids",
        output_base = data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled"
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        "plink2 --bfile {params.input_base} --keep {input.samples} --make-bed --out {params.output_base} "
        "--threads {threads}"

rule download_hollfelder_data:
    output:
        data_path + url_hollfelder.split('/')[-1]
    params:
        data = data_path,
        url = url_hollfelder
    shell:
        "wget -q -P {params.data} {params.url}"

rule extract_hollfelder_data:
    input:
        data_path + url_hollfelder.split('/')[-1]
    output:
        multiext(data_path + "Sudanfiles/Sudan",".bed",".bim",".fam")
    params:
        data_path = data_path
    shell:
        "tar xvzf {input} -C {params.data_path}"

rule select_samples_hollfelder:
    input:
        fam=data_path + "Sudanfiles/Sudan.fam"
    output:
        temp("samples_hollfelder.txt")
    params:
        exclude= "|".join(populations_to_exclude_hollfelder),
        n_samples=samples_downsampling
    shell:
        "cut -d ' ' -f1 {input.fam} | sort | uniq | xargs -I {{}} bash -c \"grep {{}} {input.fam} | cut -d ' ' -f1,2 | "
        "shuf -n {params.n_samples}\" bash | grep -v -E \"{params.exclude}\" > {output}"

rule extract_samples_hollfelder:
    input:
        to_keep = "samples_hollfelder.txt",
        plink=multiext(data_path + "Sudanfiles/Sudan", ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + "Sudanfiles/Sudan_downsampled",".bed",".bim",".fam"))
    params:
        input_base = data_path + "Sudanfiles/Sudan",
        output_base = data_path+ "Sudanfiles/Sudan_downsampled"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --keep {input.to_keep} --make-bed --set-all-var-ids @:# "
        "--out {params.output_base} --threads {threads}"

rule get_duplicate_var_ids_hollfelder:
    input:
        data_path + "Sudanfiles/Sudan_downsampled.bim"
    output:
        temp('dup_var_ids_hollfelder.txt')
    shell:
        "cut -f2 {input} | sort | uniq -d > {output}"

rule remove_duplicate_vars_hollfelder:
    input:
        plink=multiext(data_path + "Sudanfiles/Sudan_downsampled",".bed",".bim",".fam"),
        ids='dup_var_ids_hollfelder.txt'
    output:
        multiext(data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars",".bed",".bim",".fam")
    params:
        input_base=data_path + "Sudanfiles/Sudan_downsampled",
        outputbase=data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --exclude {input.ids} --make-bed --out {params.outputbase} "
        "--threads {threads}"

rule download_arauna:
    output:
        multiext(data_path + 'arauna_et_al_2017', ".bed", ".bim", '.fam')
    params:
        url_bed = url_arauna_bed,
        url_bim = url_arauna_bim,
        url_fam = url_arauna_fam,
        name_bed = url_arauna_bed.split('/')[-1],
        name_bim = url_arauna_bim.split('/')[-1],
        name_fam = url_arauna_fam.split('/')[-1],
        data_path = data_path
    shell:
        "wget -q {params.url_bed}; mv {params.name_bed} {params.data_path}arauna_et_al_2017.bed; "
        "wget -q {params.url_bim}; mv {params.name_bim} {params.data_path}arauna_et_al_2017.bim; "
        "wget -q {params.url_fam}; mv {params.name_fam} {params.data_path}arauna_et_al_2017.fam;"

# only populations with geo data available and downsample
rule select_samples_arauna:
    input:
        fam=data_path + 'arauna_et_al_2017.fam',
        coords=data_path + geo_coords_arauna
    output:
        temp("samples_to_keep_arauna.txt")
    params:
        n_samples = samples_downsampling
    shell:
        "cut -d ' ' -f1 {input.coords} | tail -n+2 | xargs -I {{}} bash -c \"grep -w {{}} {input.fam} | "
        "shuf -n{params.n_samples} | cut -d ' ' -f1,2\" > {output}"

rule extract_samples_arauna:
    input:
        to_keep="samples_to_keep_arauna.txt",
        plink=multiext(data_path + 'arauna_et_al_2017',".bim",".bed",".fam")
    output:
        multiext(data_path + 'arauna_et_al_2017_downsampled',".bim",".bed",".fam")
    params:
        input_base=data_path + 'arauna_et_al_2017',
        output_base=data_path + 'arauna_et_al_2017_downsampled'
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --keep {input.to_keep} --make-bed --set-all-var-ids @:# "
        "--out {params.output_base} --threads {threads}"

use rule select_samples_arauna as select_samples_fortes_lima with:
    input:
        fam=data_path + plink_prefix_fortes_lima + ".fam",
        coords= data_path + geo_coords_fortes_lima
    output:
        temp("samples_to_keep_fortes_lima.txt")
    params:
        n_samples = samples_downsampling

use rule extract_samples_arauna as extract_samples_fortes_lima with:
    input:
        to_keep="samples_to_keep_fortes_lima.txt",
        plink=multiext(data_path + plink_prefix_fortes_lima,".bim",".bed",".fam")
    output:
        multiext(data_path + plink_prefix_fortes_lima + '_downsampled',".bim",".bed",".fam")
    params:
        input_base=data_path + plink_prefix_fortes_lima,
        output_base=data_path + plink_prefix_fortes_lima + '_downsampled'
    threads: 32

rule get_variants_crawford_and_scheinfeldt:
    input:
        data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled.bim"
    output:
        temp("snps_crawford_and_scheinfeldt.txt")
    shell:
        "cut -f2 {input} | sort > {output}"

use rule get_variants_crawford_and_scheinfeldt as get_variants_sgdp with:
    input:
        data_path + "SGDP/downsampled_updated_ids.bim"
    output:
        temp("snps_sgdp.txt")

rule intersect_variants_sgdp_crawford_and_scheinfeldt:
    input:
        "snps_crawford_and_scheinfeldt.txt",
        "snps_sgdp.txt"
    output:
        temp("intersection_snps_sgdp_crawford_and_scheinfeldt.txt")
    shell:
        "comm -12 {input} > {output}"

rule extract_shared_snps_crawford_and_scheinfeldt:
    input:
        bed=data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled.bed",
        snps="intersection_snps_sgdp_crawford_and_scheinfeldt.txt"
    output:
        temp(multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled_intersected_snps",
            ".bed", ".bim", ".fam"))
    params:
        input_base = data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled",
        output_base = data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled_intersected_snps"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --extract {input.snps} --make-bed --out {params.output_base}"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_sgdp with:
    input:
        bed = data_path + "SGDP/downsampled_updated_ids.bed",
        snps = "intersection_snps_sgdp_crawford_and_scheinfeldt.txt"
    output:
        temp(multiext(data_path + "SGDP/downsampled_updated_ids_intersected_snps",".bed",".bim",".fam"))
    params:
        input_base=data_path + "SGDP/downsampled_updated_ids",
        output_base=data_path + "SGDP/downsampled_updated_ids_intersected_snps"

rule merge_crawford_and_scheinfeldt_and_sgdp_data:
    input:
        multiext(data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled_intersected_snps",
            ".bed", ".bim", ".fam"),
        multiext(data_path + "SGDP/downsampled_updated_ids_intersected_snps", ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + "crawford_and_scheinfeldt_sgdp_merged", ".bim", ".bed", ".fam"))
    params:
        base_a = data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled_intersected_snps",
        base_b = data_path + "SGDP/downsampled_updated_ids_intersected_snps",
        output_base = data_path + "crawford_and_scheinfeldt_sgdp_merged",
        maf = maf,
        geno = geno
    conda:
        "envs/plink.yaml"
    threads: 32
    retries: 2
    shell:
        "if [[ -f {params.output_base}_trial1-merge.missnp ]]; " # second trial --> flip SNPs once
        "then "
            "plink --bfile {params.base_b} --flip {params.output_base}_trial1-merge.missnp --make-bed --out "
            "{params.base_b}_tmp --threads {threads}; " # flip once
            "rm {params.output_base}_trial1-merge.missnp; "
            "plink --bfile {params.base_a} --bmerge {params.base_b}_tmp --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base}_trial2 --threads {threads}; "
            "mv {params.output_base}_trial2.bed {params.output_base}.bed; "
            "mv {params.output_base}_trial2.bim {params.output_base}.bim; "
            "mv {params.output_base}_trial2.fam {params.output_base}.fam; "
        "elif [[ -f {params.output_base}_trial2-merge.missnp ]]; " # last trial --> excluding SNPS that caused problems
        "then "
            "plink --bfile {params.base_a} --exclude {params.output_base}_trial2-merge.missnp --make-bed "
            "--out {params.base_a}_tmp --threads {threads}; "
            "plink --bfile {params.base_b}_tmp --exclude {params.output_base}_trial2-merge.missnp --make-bed "
            "--out {params.base_b}_tmp1 --threads {threads}; "
            "plink --bfile {params.base_a}_tmp --bmerge {params.base_b}_tmp1 --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base} --threads {threads}; "
            "rm {params.base_a}_tmp.*; "
            "rm {params.base_b}_tmp.*; "
            "rm {params.base_b}_tmp1.*; "
            "rm {params.output_base}_trial2-merge.missnp; "
        "else " # first pass
            "plink --bfile {params.base_a} --bmerge {params.base_b} --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base}_trial1 --threads {threads}; "
            "mv {params.output_base}_trial1.bed {params.output_base}.bed; "
            "mv {params.output_base}_trial1.bim {params.output_base}.bim; "
            "mv {params.output_base}_trial1.fam {params.output_base}.fam; "
        "fi"

rule get_at_cg_snps_merged:
    input:
        data_path + "crawford_and_scheinfeldt_sgdp_merged.bim"
    output:
        temp("at_gc_snps_merged_crawford_and_scheinfeldt_sgdp.txt")
    shell:
        "awk -F '\\t' '{{if ($5 == \"A\" && $6 != \"T\" || $5 != \"A\") print $0}}' {input} | "
        "awk -F '\\t' '{{if ($5 == \"T\" && $6 != \"A\" || $5 != \"T\") print $0}}' | "
        "awk -F '\\t' '{{if ($5 == \"C\" && $6 != \"G\" || $5 != \"C\") print $0}}' | "
        "awk -F '\\t' '{{if ($5 == \"G\" && $6 != \"C\" || $5 != \"G\") print $2}}' | sort > {output}"

rule filter_at_cg_snps_merged:
    input:
        "at_gc_snps_merged_crawford_and_scheinfeldt_sgdp.txt",
        multiext(data_path+ "crawford_and_scheinfeldt_sgdp_merged", ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg", ".bed", ".bim", ".fam"))
    params:
        input_base = data_path+ "crawford_and_scheinfeldt_sgdp_merged",
        output_base = data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg"
    threads: 32
    conda:
        "envs/plink.yaml"
    shell:
        "plink2 --bfile {params.input_base} --extract {input[0]} --make-bed --out {params.output_base} "
        "--threads {threads}"

rule get_positions_with_missing_ref_alt_allele_hollfelder:
    input:
        data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars.bim"
    output:
        temp('snps_missing_ref_alt_hollfelder.txt')
    shell:
        "awk  -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($5 == \".\" || $6 == \".\") print $2}}' {input} > {output}"

rule get_snps_hollfelder:
    input:
        bim=data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars.bim",
        exclude='snps_missing_ref_alt_hollfelder.txt'
    output:
        temp("snps_hollfelder.txt")
    shell:
        "grep -v -w -f {input.exclude} {input.bim} | cut -f2 | sort > {output}"

rule intersect_snps_sgdp_crawford_and_scheinfeldt_and_hollfelder:
    input:
        "snps_hollfelder.txt",
        "at_gc_snps_merged_crawford_and_scheinfeldt_sgdp.txt"
    output:
        temp("intersection_snps_crawford_and_scheinfeldt_sgdp_hollfelder.txt")
    shell:
        "comm -12 {input} > {output}"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_crawford_and_scheinfeldt_sgdp with:
    input:
        plink = multiext(data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg", ".bed", ".bim", ".fam"),
        snps = "intersection_snps_crawford_and_scheinfeldt_sgdp_hollfelder.txt"
    output:
        temp(multiext(data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg_intersected_snps_hollfelder",
                      ".bed",".bim",".fam"))
    params:
        input_base=data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg",
        output_base=data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg_intersected_snps_hollfelder"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_hollfelder with:
    input:
        bed = data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars.bed",
        snps = "intersection_snps_crawford_and_scheinfeldt_sgdp_hollfelder.txt"
    output:
        temp(multiext(data_path +
                      "Sudanfiles/Sudan_downsampled_rm_dup_vars_intersected_snps_crawford_and_scheinfeldt_sgdp",
                      ".bed",".bim",".fam"))
    params:
        input_base=data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars",
        output_base=data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars_intersected_snps_crawford_and_scheinfeldt_sgdp"

use rule merge_crawford_and_scheinfeldt_and_sgdp_data as merge_crawford_and_scheinfeldt_sgdp_and_hollfelder_data with:
    input:
        multiext(data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg_intersected_snps_hollfelder",
                 ".bed", ".bim", ".fam"),
        multiext(data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars_intersected_snps_crawford_and_scheinfeldt_sgdp",
                 ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged", ".bim", ".bed", ".fam"))
    params:
        base_a = data_path+ "crawford_and_scheinfeldt_sgdp_merged_filtered_at_cg_intersected_snps_hollfelder",
        base_b = data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars_intersected_snps_crawford_and_scheinfeldt_sgdp",
        output_base = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged",
        maf= maf,
        geno= geno
    threads: 32

rule get_snps_arauna:
    input:
        data_path + "arauna_et_al_2017_downsampled.bim"
    output:
        temp('snps_arauna.txt')
    shell:
        'cut -f2 {input} | sort > {output}'

rule get_snps_crawford_and_scheinfeldt_sgdp_hollfelder_merged:
    input:
        data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged.bim"
    output:
        temp('snps_crawford_and_scheinfeldt_sgdp_hollfelder_merged.txt')
    shell:
        'cut -f2 {input} | sort > {output}'

rule intersect_snps_arauna_and_crawford_and_scheinfeldt_sgdp_hollfelder_merged:
    input:
        'snps_arauna.txt',
        'snps_crawford_and_scheinfeldt_sgdp_hollfelder_merged.txt'
    output:
        temp('intersect_snps_arauna_and_crawford_and_scheinfeldt_sgdp_hollfelder_merged.txt')
    shell:
        "comm -12 {input} > {output}"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_arauna_and_crawford_and_scheinfeldt_sgdp_hollfelder with:
    input:
        bed = data_path + "arauna_et_al_2017_downsampled.bed",
        snps = 'intersect_snps_arauna_and_crawford_and_scheinfeldt_sgdp_hollfelder_merged.txt'
    output:
        temp(multiext(data_path +
                      "arauna_et_al_2017_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder",
                      ".bed",".bim",".fam"))
    params:
        input_base=data_path + "arauna_et_al_2017_downsampled",
        output_base=data_path + "arauna_et_al_2017_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_crawford_and_scheinfeldt_sgdp_hollfelder_and_arauna with:
    input:
        plink = multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged", ".bim", ".bed", ".fam"),
        snps = 'intersect_snps_arauna_and_crawford_and_scheinfeldt_sgdp_hollfelder_merged.txt'
    output:
        temp(multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged_intersected_snps_arauna",
             ".bed",".bim",".fam"))
    params:
        input_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged",
        output_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged_intersected_snps_arauna"

use rule merge_crawford_and_scheinfeldt_and_sgdp_data as merge_crawford_and_scheinfeldt_sgdp_hollfelder_and_arauna_data with:
    input:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged_intersected_snps_arauna", ".bed", ".bim", ".fam"),
        multiext(data_path + "arauna_et_al_2017_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder",
                 ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged", ".bim", ".bed", ".fam"))
    params:
        base_a = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_merged_intersected_snps_arauna",
        base_b = data_path + "arauna_et_al_2017_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder",
        output_base = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged",
        maf= maf,
        geno= geno
    threads: 32

rule get_snps_fortes_lima:
    input:
        data_path + plink_prefix_fortes_lima + "_downsampled.bim"
    output:
        temp("snps_fortes_lima.txt")
    shell:
        "cut -f2 {input} | sort > {output}"

rule get_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna:
    input:
        data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged.bim"
    output:
        temp("snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged.txt")
    shell:
        "cut -f2 {input} | sort > {output}"

rule intersect_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima:
    input:
        "snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged.txt",
        "snps_fortes_lima.txt"
    output:
        temp("intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima.txt")
    shell:
        "comm -12 {input} > {output}"

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_fortes_lima_and_crawford_and_scheinfeldt_sgdp_hollfelder_arauna with:
    input:
        bed = data_path + plink_prefix_fortes_lima + "_downsampled.bed",
        snps = 'intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima.txt'
    output:
        temp(multiext(data_path + plink_prefix_fortes_lima +
                      "_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna",
            ".bed",".bim",".fam"))
    params:
        input_base=data_path + plink_prefix_fortes_lima + "_downsampled",
        output_base=data_path + plink_prefix_fortes_lima + "_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna"
    threads: 32

use rule extract_shared_snps_crawford_and_scheinfeldt as extract_shared_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima with:
    input:
        plink = multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged", ".bim", ".bed", ".fam"),
        snps = 'intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima.txt'
    output:
        temp(multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged_intersected_snps_fortes_lima",
                      ".bed",".bim",".fam"))
    params:
        input_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged",
        output_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged_intersected_snps_fortes_lima"
    threads: 32

use rule merge_crawford_and_scheinfeldt_and_sgdp_data as merge_crawford_and_scheinfeldt_sgdp_hollfelder_arauna_and_fortes_lima_data with:
    input:
        multiext(data_path + plink_prefix_fortes_lima +
                 "_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna",
                 ".bed", ".bim", ".fam"),
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged_intersected_snps_fortes_lima",
                 ".bed", ".bim", ".fam")
    output:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
                 ".bim", ".bed", ".fam")
    params:
        base_a = data_path + plink_prefix_fortes_lima + "_downsampled_intersected_snps_crawford_and_scheinfeldt_sgdp_hollfelder_arauna",
        base_b = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_merged_intersected_snps_fortes_lima",
        output_base = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
        maf= maf,
        geno= geno
    threads: 32

rule find_snps_in_ld:
    input:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
                 ".bed", ".bim", ".fam")
    output:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
                 ".prune.in", ".prune.out",)
    params:
        input_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
        max_r2 = max_r2
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        "plink2 --bfile {params.input_base} --indep-pairwise 50 1 {params.max_r2} --out {params.input_base} "
        "--threads {threads}"

rule ld_prune:
    input:
        extract=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged.prune.in",
        bed=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged.bed"
    output:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned",
                 ".bed", ".bim", ".fam")
    params:
        input_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged",
        output_base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned"
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        "plink2 --bfile {params.input_base} --extract {input.extract} --make-bed --out {params.output_base} "
        "--threads {threads}"

rule perform_pca:
    input:
        data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.bed"
    output:
        multiext(data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned",
                 ".eigenvec", ".eigenval")
    params:
        base=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned"
    conda:
        "envs/plink.yaml"
    threads: 32
    shell:
        "plink2 --bfile {params.base} --pca 20 --out {params.base} --threads {threads}"

rule plot_pca:
    input:
        eigenvec = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.eigenvec",
        coords = data_path + "individual_geo_coords.tab"
    output:
        results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_pca.png"
    params:
        output_dir = results_path
    shell:
        "scripts/plot_pca.py --eigenvec {input.eigenvec} --plot_dir {params.output_dir} --coords {input.coords}"


rule get_individual_coordinates:
    input:
        meta_sgdp = data_path + "SGDP/" + url_metadata_sgdp.split('/')[-1],
        meta_crawford_and_scheinfeldt = data_path + geo_coords_crawford_and_scheinfeldt,
        meta_hollfelder = data_path + geo_coords_hollfelder,
        meta_arauna = data_path + geo_coords_arauna,
        meta_fortes_lima = data_path + geo_coords_fortes_lima,
        sgdp_fam = data_path + "SGDP/downsampled_updated_ids.fam",
        crawford_and_scheinfeldt_fam = data_path + plink_prefix_crawford_and_scheinfeldt + "_downsampled.fam",
        hollfelder_fam = data_path + "Sudanfiles/Sudan_downsampled_rm_dup_vars.fam",
        arauna_fam = data_path + "arauna_et_al_2017_downsampled.fam",
        fortes_lima_fam = data_path + plink_prefix_fortes_lima + "_downsampled.fam"
    output:
        data_path + "individual_geo_coords.tab"
    shell:
        "cut -f2 {input.sgdp_fam} | xargs -I {{}} bash -c \"grep {{}} {input.meta_sgdp} | cut -f6,7\" | "
        "paste <(cut -f1,2 {input.sgdp_fam}) - > {output}; "
        "cut -f2 {input.crawford_and_scheinfeldt_fam} | "
        "xargs -I {{}} bash -c \"grep {{}} {input.meta_crawford_and_scheinfeldt}\" | cut -f4,5 | "
        "paste <(cut -f1,2 {input.crawford_and_scheinfeldt_fam}) - >> {output}; "
        "cut -f1 {input.hollfelder_fam} | xargs -I {{}} bash -c \"grep {{}} {input.meta_hollfelder}\" | cut -f6,7 | "
        "paste <(cut -f1,2 {input.hollfelder_fam}) - >> {output}; "
        "cut -f1 {input.arauna_fam} | xargs -I {{}} bash -c \"grep -w {{}} {input.meta_arauna}\" | cut -d ' ' -f2,3 | "
        "paste <(cut -f1,2 {input.arauna_fam}) - | sed -e 's/\\s/\\t/g' >> {output}; "
        "cut -f1 {input.fortes_lima_fam} | xargs -I {{}} bash -c \"grep -w {{}} {input.meta_fortes_lima}\" | "
        "cut -d ' ' -f2,3 | paste <(cut -f1,2 {input.fortes_lima_fam}) - | sed -e 's/\\s/\\t/g' >> {output}"

rule run_admixture:
    input:
        data_path +  "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.bed"
    output:
        multiext(results_path +  "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.{K}",
                 '.Q', '.P', '.log')
    params:
        data = data_path,
        results = results_path,
        prefix = "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned"
    threads: 32
    conda:
        "envs/admixture.yaml"
    shell:
        "admixture --cv {input} {wildcards.K} -j{threads} > {params.prefix}.{wildcards.K}.log; "
        "mkdir -p {params.results}; "
        "mv {params.prefix}.{wildcards.K}.* {params.results}"

rule zip_q_files:
    input:
        expand(results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.{K}.Q",K=K)
    output:
        results_path + "qfiles.zip"
    params:
        files = expand("crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.{K}.Q",K=K),
        results_path = results_path,
        zip = "qfiles.zip"
    shell:
        "cd {params.results_path}; zip {params.zip} {params.files}"

rule make_population_file_clumpak:
    input:
        data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.fam"
    output:
        temp('populations.txt')
    shell:
        "cut -f1 {input} > {output}"

rule align_admixture_results_across_k:
    input:
        clumpak_dir = url_clumpak.split('/')[-1].split('.zip')[0],
        qfiles = results_path + "qfiles.zip",
        pops = 'populations.txt'
    output:
        directory('results/clumpak')
    conda:
        "envs/clumpak.yaml"
    shell:
        "cd {input.clumpak_dir}; cpanm List::Permutor; "
        "perl distructForManyKs.pl --id 1 --dir ../{output} --file ../{input.qfiles} "
        "--inputtype admixture --indtopop ../{input.pops}"

rule visualize_all_admixture_results:
    input:
        clumpak = 'results/clumpak',
        fam = data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.fam",
        coords= data_path + "individual_geo_coords.tab"
    output:
        all=results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_plots.pdf",
        best_k=results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_plot_best_K.pdf",
        kriging=results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_admixture_kriging_best_K.pdf"
    params:
        results_path = results_path,
        output_prefix = results_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned"
    conda:
        "envs/admixture.yaml"
    shell:
        "scripts/visualize_admixture.py -i {input.clumpak}/aligned.files/ -f {input.fam} -l {params.results_path} "
        "-c {input.coords} -o {params.output_prefix}"

rule run_feems:
    input:
        pop_coords = data_path + "individual_geo_coords.tab",
        bed=data_path + "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned.bed"
    output:
        results_path +  "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned_feems_plot.pdf"
    params:
        data_path = data_path,
        prefix = "crawford_and_scheinfeldt_sgdp_hollfelder_arauna_fortes_lima_merged_ld_pruned",
        output_path = results_path
    threads: 1
    conda:
        "envs/feems_e.yaml"
    shell:
        "scripts/run_feems.py -d {params.data_path} -p {params.prefix} -c {input.pop_coords} -o {params.output_path}"


