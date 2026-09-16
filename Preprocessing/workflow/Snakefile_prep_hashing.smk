"""Helpers for the HASHED (multiplexed nuclei) samples.

Mirrors Snakefile_prep_RNA.smk, but each 10x channel carries two libraries --
a Gene Expression library and an Antibody Capture (TotalSeq-B hashtag) library
-- which `cellranger count` must process together via a libraries.csv.
"""
import glob
import os
import re
import sys

samples_list = list(config['samples'].keys())
print(samples_list)


def parse_fastq(path):
    """Split a 10x FASTQ filename into its Illumina fields."""
    result = {}
    fastq = os.path.basename(path)
    result['number'] = re.findall('_S[0-9]+_', fastq)[0].strip('_')
    result['lane'] = re.findall('_L[0-9]+_', fastq)[0].strip('_')
    result['read'] = re.findall('_[RI][0-9]+_', fastq)[0].strip('_')
    result['id'] = re.split('_S[0-9]+_', fastq)[0].strip('_')
    result['suffix'] = re.split('_[RI][0-9]+_', fastq)[1].strip('_')
    return result


def list_fastq(fastq_folder):
    """All FASTQs under a library folder, searched recursively."""
    files = sorted(glob.glob(os.path.join(fastq_folder, '**', '*.fastq.gz'),
                             recursive=True))
    if not files:
        sys.stderr.write('*** Error: Found 0 fastq files in {}\n'.format(fastq_folder))
        raise Exception('No files found in fastq folder\n')
    return files


def get_fastq_for_sample(wildcards):
    """Both libraries of a channel are inputs to its cellranger run."""
    entry = config['samples'][wildcards.sample]
    return (list_fastq(entry['fastq_path_RNA'])
            + list_fastq(entry['fastq_path_feature']))


def fastq_sample_prefixes(fastq_folder):
    """The distinct `--sample` prefixes present in a library folder.

    Cell Ranger matches FASTQs by this prefix, which is set at demultiplexing
    and is generally NOT the same string for the GEX and the feature library of
    one channel -- hence it is read off the filenames rather than assumed.
    """
    prefixes = sorted({parse_fastq(f)['id'] for f in list_fastq(fastq_folder)})
    if not prefixes:
        raise Exception('Could not parse any sample prefix in {}\n'.format(fastq_folder))
    return prefixes


def write_libraries_csv(sample, path):
    """Emit the libraries.csv consumed by `cellranger count --libraries`."""
    entry = config['samples'][sample]
    rows = []
    for key, library_type in (('fastq_path_RNA', 'Gene Expression'),
                              ('fastq_path_feature', 'Antibody Capture')):
        folder = os.path.abspath(entry[key])
        for prefix in fastq_sample_prefixes(folder):
            rows.append((folder, prefix, library_type))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write('fastqs,sample,library_type\n')
        for row in rows:
            handle.write('{0},{1},{2}\n'.format(*row))
    return rows
