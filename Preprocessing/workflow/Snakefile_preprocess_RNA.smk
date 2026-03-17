include: 'Snakefile_prep_RNA.smk'

rule all_preprocess:
    input:
        cellbender_seurat_h5=[
            'Humanized_mice_OD/{sample}/RNA_AAAAGGGG/epochs_300/cellbender_output_seurat.h5'.format(sample=sample) for sample in samples_list]
       

rule run_cellranger_RNA:
    input:
        lambda wildcards: get_fastq_for_cellranger_rna(config['samples'][wildcards.sample]['fastq_path_RNA']+ '/**/*{lane}*R[12]*.fastq.gz',sample=wildcards.sample)
    output:
        bam_RNA='Humanized_mice_OD/{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam'
    params:
        cellranger_software=config['general']['cellranger_software_RNA'],
        cellranger_ref=config['general']['cellranger_ref_RNA'],
        fastq_folder=lambda wildcards: config['samples'][wildcards.sample]['fastq_path_RNA'],
        mem=64
    threads: 20
    resources:
        mem_mb = 64000
    shell:
        'mkdir -p {wildcards.sample}/RNA_AAAAGGGG; '
        'cd {wildcards.sample}/RNA_AAAAGGGG/; '
        'rm -rf cellranger/; '
        '{params.cellranger_software} count --id cellranger --transcriptome {params.cellranger_ref}  --fastqs {params.fastq_folder} --localcores={threads} --localmem={params.mem} --create-bam=true'


rule run_cellbender:
    input:
        h5='Humanized_mice_OD/{sample}/RNA_AAAAGGGG/cellranger/outs/raw_feature_bc_matrix.h5'
    output:
        cellbender_h5='Humanized_mice_OD/{sample}/RNA_AAAAGGGG/epochs_300/cellbender_output.h5'
    params:
        fpr=config['general']['cellbender_fpr'],
        cuda_flag=lambda wildcards: '--cuda' if config['general']['cellbender_cuda'] else '',
        epochs=config['general']['cellbender_epochs']
    threads: 1
    resources:
        mem_mb = 16000,
        gpu=1
    conda:
        '/home/mattia/miniconda3/envs/cellbender.yml'
    shell:
        'cellbender remove-background '
        '--input {input.h5} '
        '--output {output.cellbender_h5} '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.cuda_flag}'


rule ptrepack_seurat:
    input:
        cellbender_h5='Humanized_mice_OD/{sample}/RNA_AAAAGGGG/epochs_300/cellbender_output.h5'
    output:
        seurat_h5='Humanized_mice_OD/{sample}/RNA_AAAAGGGG/epochs_300/cellbender_output_seurat.h5'
    threads: 1
    resources:
        mem_mb = 16000
    conda:
        '/home/mattia/miniconda3/envs/cellbender.yml'
    shell:
        'ptrepack --complevel 5 {input.cellbender_h5}:/matrix {output.seurat_h5}:/matrix'
