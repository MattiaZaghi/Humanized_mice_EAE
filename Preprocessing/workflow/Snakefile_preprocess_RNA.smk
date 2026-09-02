include: 'Snakefile_prep_RNA.smk'

rule all_preprocess:
    input:
        cellbender_seurat_h5=[
            '{sample}/RNA_AAAAGGGG/epochs/cellbender_output_seurat.h5'.format(sample=sample) for sample in samples_list]
       

rule run_cellranger_RNA:
    input:
        lambda wildcards: get_fastq_for_cellranger_rna(config['samples'][wildcards.sample]['fastq_path_RNA']+ '/**/*{lane}*R[12]*.fastq.gz',sample=wildcards.sample)
    output:
        bam_RNA='{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam',
        h5='{sample}/RNA_AAAAGGGG/cellranger/outs/raw_feature_bc_matrix.h5'
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
        h5='{sample}/RNA_AAAAGGGG/cellranger/outs/raw_feature_bc_matrix.h5'
    output:
        cellbender_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output.h5'
    params:
        fpr=config['general']['cellbender_fpr'],
        cuda_flag=lambda wildcards: '--cuda' if config['general']['cellbender_cuda'] else '',
        epochs=config['general']['cellbender_epochs']
    threads: 1
    resources:
        mem_mb = 16000,
        gpu=1
    conda:
        '/home/mattia/miniconda3_n/envs/cellbender.yml'
    shell:
        # run via the patched wrapper: strips pyro's un-pickleable .unconstrained
        # weakref off model params before CellBender's torch.save checkpoint, which
        # otherwise crashes with "cannot pickle 'weakref.ReferenceType' object" and
        # aborts the run (see workflow/scripts/run_cellbender_patched.py).
        'python /home/mattia/Humanized_mice_EAE/Preprocessing/workflow/scripts/run_cellbender_patched.py remove-background '
        '--input {input.h5} '
        '--output {output.cellbender_h5} '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.cuda_flag}'


rule ptrepack_seurat:
    input:
        cellbender_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output.h5'
    output:
        seurat_h5='{sample}/RNA_AAAAGGGG/epochs/cellbender_output_seurat.h5'
    threads: 1
    resources:
        mem_mb = 16000
    conda:
        '/home/mattia/miniconda3_n/envs/cellbender.yml'
    shell:
        # ptrepack aborts with "node names ... duplicated in destination" if the
        # output already exists (e.g. a partial file left by a killed/retried job),
        # causing a permanent failure loop on reruns. Remove any stale destination
        # first and pass --overwrite-nodes so the step is idempotent.
        # Also call the conda env's ptrepack explicitly ($CONDA_PREFIX/bin): a bare
        # `ptrepack` resolves to ~/.local/bin/ptrepack (system python, no `tables`)
        # -> ModuleNotFoundError: No module named 'tables'. PYTHONNOUSERSITE=1 keeps
        # the user site-packages out of the env python.
        'rm -f {output.seurat_h5}; '
        'PYTHONNOUSERSITE=1 "$CONDA_PREFIX/bin/ptrepack" --overwrite-nodes --complevel 5 {input.cellbender_h5}:/matrix {output.seurat_h5}:/matrix'
