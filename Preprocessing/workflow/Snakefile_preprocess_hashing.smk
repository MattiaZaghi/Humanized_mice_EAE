include: 'Snakefile_prep_hashing.smk'

rule all_preprocess_hashing:
    input:
        cellbender_seurat_h5=[
            '{sample}/RNA_AAAAGGGG/epochs/cellbender_output_seurat.h5'.format(sample=sample)
            for sample in samples_list]


rule make_libraries_csv:
    """One CSV per channel listing its GEX and Antibody Capture libraries."""
    output:
        libraries='{sample}/RNA_AAAAGGGG/libraries.csv'
    run:
        rows = write_libraries_csv(wildcards.sample, output.libraries)
        for row in rows:
            print('{0}: {2} <- {1}'.format(wildcards.sample, row[1], row[2]))


rule run_cellranger_hashing:
    input:
        fastqs=get_fastq_for_sample,
        libraries='{sample}/RNA_AAAAGGGG/libraries.csv'
    output:
        bam_RNA='{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam',
        h5='{sample}/RNA_AAAAGGGG/cellranger/outs/raw_feature_bc_matrix.h5'
    params:
        cellranger_software=config['general']['cellranger_software_RNA'],
        cellranger_ref=config['general']['cellranger_ref_RNA'],
        feature_ref=config['general']['feature_reference'],
        chemistry=config['general'].get('chemistry', 'auto'),
        mem=64
    threads: 20
    resources:
        mem_mb=64000
    shell:
        # --libraries replaces --fastqs: it points cellranger at BOTH the
        # transcriptome and the hashtag library, and --feature-ref declares the
        # four TotalSeq-B hashtag barcodes. The raw_feature_bc_matrix.h5 then
        # contains a "Gene Expression" and an "Antibody Capture" feature block.
        'mkdir -p {wildcards.sample}/RNA_AAAAGGGG; '
        'cd {wildcards.sample}/RNA_AAAAGGGG/; '
        'rm -rf cellranger/; '
        '{params.cellranger_software} count --id cellranger '
        '--transcriptome {params.cellranger_ref} '
        '--libraries $(realpath ../../{input.libraries}) '
        '--feature-ref {params.feature_ref} '
        '--chemistry {params.chemistry} '
        '--localcores={threads} --localmem={params.mem} --create-bam=true'


rule run_cellbender_hashing:
    input:
        h5='{sample}/RNA_AAAAGGGG/cellranger/outs/raw_feature_bc_matrix.h5'
    output:
        cellbender_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output.h5'
    params:
        fpr=config['general']['cellbender_fpr'],
        cuda_flag=lambda wildcards: '--cuda' if config['general']['cellbender_cuda'] else '',
        epochs=config['general']['cellbender_epochs']
    threads: 1
    resources:
        mem_mb=16000,
        gpu=1
    conda:
        '/home/mattia/miniconda3_n/envs/cellbender.yml'
    shell:
        # Same patched wrapper as the unhashed samples (see
        # workflow/scripts/run_cellbender_patched.py for the weakref fix).
        # CellBender denoises the Antibody Capture block alongside the genes, so
        # the hashtag counts that reach Seurat are already ambient-corrected.
        'python /home/mattia/Humanized_mice_EAE/Preprocessing/workflow/scripts/run_cellbender_patched.py remove-background '
        '--input {input.h5} '
        '--output {output.cellbender_h5} '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.cuda_flag}'


rule ptrepack_seurat_hashing:
    input:
        cellbender_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output.h5'
    output:
        seurat_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output_seurat.h5'
    threads: 1
    resources:
        mem_mb=16000
    conda:
        '/home/mattia/miniconda3_n/envs/cellbender.yml'
    shell:
        'rm -f {output.seurat_h5}; '
        'PYTHONNOUSERSITE=1 "$CONDA_PREFIX/bin/ptrepack" --overwrite-nodes --complevel 5 '
        '{input.cellbender_h5}:/matrix {output.seurat_h5}:/matrix'
