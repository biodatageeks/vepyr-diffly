#!/usr/bin/env python3
"""Run full VEP benchmark on HG002 4.1M records with timing and auto-cleanup."""

import os
from pathlib import Path
import shlex
import subprocess
import time

ROOT = Path('/Users/lukaszjezapkowicz/workspace/data_vepyr')
RUN = ROOT / 'runs/chapter5-m5-20260924T061203Z'
OUT = RUN / 'benchmarks/vep'
LOG = OUT / 'commands.log'
IMAGE = 'ensemblorg/ensembl-vep:release_116.0'

output_vcf = OUT / 'HG002_vep_forknone.vcf'
timing_file = OUT / 'full_forknone.time.txt'
stdout_file = OUT / 'full_forknone.stdout.txt'
stderr_file = OUT / 'full_forknone.stderr.txt'

vep_args = [
    'vep', '--dir', '/opt/vep/.vep', '--cache', '--merged', '--offline',
    '--assembly', 'GRCh38',
    '--input_file', '/input/HG002_normalized.vcf.gz',
    '--output_file', '/output/HG002_vep_forknone.vcf',
    '--vcf', '--force_overwrite', '--no_stats', '--everything', '--hgvs',
    '--fasta', '/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
]

docker_cmd = [
    'docker', 'run', '--rm', '--user', '501:20', '--env', 'HOME=/tmp',
    '-v', f'{ROOT}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro',
    '-v', f'{ROOT}/input:/input:ro',
    '-v', f'{OUT}:/output',
    IMAGE,
    *vep_args,
]

full_command = ['/opt/homebrew/bin/rtk', '/usr/bin/time', '-p', '-o', str(timing_file), *docker_cmd]

print('Logging command to:', LOG, flush=True)
with LOG.open('a') as fh:
    fh.write('full_forknone: ' + shlex.join(full_command) + '\n')

print('Starting full VEP benchmark on 4,096,123 variants...', flush=True)
started_at = time.time()
with stdout_file.open('w') as out_fh, stderr_file.open('w') as err_fh:
    res = subprocess.run(full_command, stdout=out_fh, stderr=err_fh, check=False)

elapsed = time.time() - started_at
print(f'Full VEP finished with exit code {res.returncode} in {elapsed:.2f} s', flush=True)

if res.returncode != 0:
    print('Error occurred! See', stderr_file, flush=True)
    raise SystemExit(res.returncode)

if output_vcf.exists():
    print('Counting lines in output VCF...', flush=True)
    with output_vcf.open('rb') as fh:
        lines = sum(1 for _ in fh)
    print(f'Lines: {lines}', flush=True)
    if timing_file.exists():
        print(f'Timing:\n{timing_file.read_text().strip()}', flush=True)
    print(f'Removing large VCF ({output_vcf.stat().st_size / 1024**3:.2f} GB) to free disk space...', flush=True)
    output_vcf.unlink()
    print('VCF removed. Benchmark complete.', flush=True)
