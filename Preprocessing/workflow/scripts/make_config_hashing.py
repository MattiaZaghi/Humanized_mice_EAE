#!/usr/bin/env python
"""Build config_hashing.yaml from the raw hashing directory.

The hashed run stores two FASTQ folders per 10x channel -- one whose name
contains "RNA" (transcriptome) and one whose name contains "feature" (TotalSeq-B
hashtag library). This script pairs them by the channel name that remains once
the RNA/feature token is stripped, so the config never has to be typed by hand.

Usage:
    python make_config_hashing.py \
        --raw-dir /ws/gcb/MZ/Humanized_mice_OD/Raw_data/hashing \
        --output  Preprocessing/config/config_hashing.yaml
"""
import argparse
import os
import re
import sys

GENERAL = {
    'cellranger_software_RNA': '/home/mattia/cellranger-9.0.1/bin/cellranger',
    'cellranger_ref_RNA': '/date/gcb/gcb_MZ/refdata-gex-GRCh38-and-mm10-2020-A/',
    'feature_reference': ('/home/mattia/Humanized_mice_EAE/Preprocessing/config/'
                          'feature_reference_TotalSeqB_NPC_hashtags.csv'),
    'chemistry': 'auto',
    'cellbender_fpr': 0.01,
    'cellbender_cuda': True,
    'cellbender_epochs': 150,
}

RNA_TOKEN = re.compile(r'[._-]?RNA[._-]?', re.IGNORECASE)
FEATURE_TOKEN = re.compile(r'[._-]?(feature|FB|HTO|ADT)[._-]?', re.IGNORECASE)


def channel_name(folder, token):
    """Strip the library-type token to get the shared channel name."""
    return token.sub('_', folder).strip('_') or folder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    folders = sorted(f for f in os.listdir(args.raw_dir)
                     if os.path.isdir(os.path.join(args.raw_dir, f)))

    rna, feature = {}, {}
    for folder in folders:
        if FEATURE_TOKEN.search(folder):
            feature[channel_name(folder, FEATURE_TOKEN)] = folder
        elif RNA_TOKEN.search(folder):
            rna[channel_name(folder, RNA_TOKEN)] = folder
        else:
            sys.stderr.write('*** Skipping unrecognised folder: {}\n'.format(folder))

    unpaired = set(rna) ^ set(feature)
    if unpaired:
        sys.stderr.write('*** Error: channels without both libraries: {}\n'
                         .format(sorted(unpaired)))
        raise SystemExit(1)
    if not rna:
        sys.stderr.write('*** Error: no channels found in {}\n'.format(args.raw_dir))
        raise SystemExit(1)

    lines = ['samples:']
    for channel in sorted(rna):
        lines.append('  {}:'.format(channel))
        lines.append('    fastq_path_RNA: {}/'.format(
            os.path.join(args.raw_dir, rna[channel])))
        lines.append('    fastq_path_feature: {}/'.format(
            os.path.join(args.raw_dir, feature[channel])))
    lines.append('')
    lines.append('general:')
    for key, value in GENERAL.items():
        lines.append('  {}: {}'.format(key, str(value).lower()
                                       if isinstance(value, bool) else value))
    lines.append('')

    with open(args.output, 'w') as handle:
        handle.write('\n'.join(lines))
    print('Wrote {} channels to {}'.format(len(rna), args.output))
    for channel in sorted(rna):
        print('  {}: {} + {}'.format(channel, rna[channel], feature[channel]))


if __name__ == '__main__':
    main()