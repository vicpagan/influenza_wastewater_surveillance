# influenza_wastewater_surveillance

Estimates the relative proportions of influenza strains in a wastewater sequencing sample.

Reads are aligned to a few reference strains, projected onto a multiple sequence alignment (MSA) of known influenza genomes, and compared against every strain in the MSA. An expectation–maximization (EM) mixture model then estimates what proportion of the sample each strain makes up.

## How it works

`eliminate_strains` runs these steps in order:

1. **Read the MSA.** Loads every strain's aligned sequence and name.
2. **Locate the reference strains.** For each strain given with `-r`, reads `<reference-sequences-dir>/<STRAIN>.fasta`, finds that strain in the MSA, and loads its base-to-column map from the positions file (`-p`).
3. **Align reads with Bowtie2.** For each reference strain, builds a Bowtie2 index in `<bowtie2-indexes-dir>/<STRAIN>` if one doesn't exist yet, then aligns the reads with `bowtie2 --all`, writing `<working-dir>/<STRAIN>.sam`.
4. **Read the SAM files and drop unaligned reads.** A read (or read pair) is dropped only if it failed to align to every reference strain. For paired reads, both mates must be mapped to count as aligned.
5. **Build the mismatch matrix.** For each read, its best alignment (fewest mismatches to its reference) is projected onto MSA columns, and its mismatches against every strain in the MSA are counted.
6. **Group indistinguishable strains.** Strains whose mismatch counts are identical across every read cannot be told apart by this sample, so they are merged into a single `Group N`.
7. **Estimate proportions with EM.** Each read's likelihood under strain *j* is

   L<sub>ij</sub> = ε<sup>m<sub>ij</sub></sup> · (1 − ε)<sup>b<sub>i</sub> − m<sub>ij</sub></sup>

   where m<sub>ij</sub> is the read's mismatch count against strain *j*, b<sub>i</sub> is its number of aligned bases, and ε is the error rate (`-e`). EM, accelerated with SQUAREM and parallelized with OpenMP, finds the strain proportions that maximize the total likelihood. It runs from two starting points, one informed by the data and one random, and prints how far apart the two answers are as a convergence check.
8. **Optional likelihood-ratio test** (`-l`). For the top `-z` strains, refits the model without each strain and reports a log-likelihood ratio measuring how much that strain is needed.

## Requirements

- **C compiler with OpenMP**, such as GCC.
- **zlib** (`-lz`). All input files may be gzip-compressed.
- **Bowtie2**. `bowtie2` and `bowtie2-build` must be on your `PATH`.

## Building

```bash
cd eliminate_strains
make
make debug
make clean
```

## Inputs

| Input | Flag | Format |
|---|---|---|
| MSA | `-i` | FASTA. Every sequence must be on **a single line**, all sequences must be the same length, and each line must be under 30,000 characters. Sequence names must match the reference strain names exactly. |
| Positions file | `-p` | For each reference strain, a header line `>STRAIN` followed by one line of comma-separated integers: the MSA column index of each base of that strain's reference sequence, in order. The number of entries must equal the reference sequence length. |
| Reference sequences | `-g` | A directory containing `<STRAIN>.fasta` for every strain passed to `-r`. The header must be the strain name, and the sequence must be on **a single line**. |
| Bowtie2 indexes | `-b` | A directory for `<STRAIN>.1.bt2` and the other index files. Missing indexes are built automatically, but the directory must already exist. |
| Reads | `-1`/`-2` or `-0` | FASTQ by default, or FASTA with `-a`. Use `-P` with `-1` and `-2` for paired-end reads, or `-0` for single-end reads. |

## Usage

```bash
./eliminate_strains \
  -i msa.fasta \
  -p strain_positions.txt \
  -g reference_sequences \
  -b bowtie2_indexes \
  -r STRAIN_1 STRAIN_2 STRAIN_3 \
  -P \
  -1 forward_reads.fastq \
  -2 reverse_reads.fastq \
  -w working \
  -o output \
  -m 250000 \
  -t 32
```

`-r` takes any number of strain names as **separate arguments**, up to the next argument that starts with `-`. Quoting them into one argument (`-r "A B C"`) makes the program look for a single strain named `A B C`.

### Options

| Short | Long | Argument | Default | Description |
|---|---|---|---|---|
| `-i` | `--MSA-filepath` | file | required | MSA of influenza strains (FASTA) |
| `-p` | `--non-imputed-positions-filepath` | file | required | Base-to-MSA-column map for each reference strain |
| `-g` | `--reference-sequences-dir` | dir | required | Directory of `<STRAIN>.fasta` reference sequences |
| `-b` | `--bowtie2-indexes-dir` | dir | required | Directory of Bowtie2 indexes, built if missing |
| `-r` | `--reference-strains` | names | required | One or more reference strain names, space-separated |
| `-P` | `--paired` | | off | Reads are paired-end |
| `-1` | `--forward-end-filepath` | file | | Forward reads (paired-end) |
| `-2` | `--reverse-end-filepath` | file | | Reverse reads (paired-end) |
| `-0` | `--single-end-filepath` | file | | Reads (single-end) |
| `-a` | `--fasta` | | FASTQ | Reads are FASTA |
| `-o` | `--output-directory` | dir | `.` | Where output files are written |
| `-w` | `--working-directory` | dir | `.` | Where SAM files are written |
| `-m` | `--max-num-reads` | int | 100000 | Stop with an error if more aligned reads than this remain (see below) |
| `-e` | `--EM-error` | decimal | 0.005 | Error rate ε used in the EM likelihood, between 0.001 and 1 |
| `-t` | `--cores` | int | 1 | Number of threads |
| `-l` | `--llr` | | off | Run the per-strain likelihood-ratio test |
| `-z` | `--num-top-strains-llr` | int | 10 | Number of top strains to compute the LLR for |
| `-n` | `--no-read-sam` | | off | Currently has no effect |
| `-v` | `--verbose` | | off | Show Bowtie2's own output instead of running it with `--quiet` |
| `-h` | `--help` | | | Print usage |

## Outputs

All output files go to the `-o` directory and are overwritten on every run.

**`proportions.csv`** lists every strain or group with its estimated proportion, sorted from highest to lowest:

```csv
names.arg,p
"STRAIN_1",0.412
"Group 0",0.301
"STRAIN_2",0.287
```

**`unidentifiable_strains.txt`** lists the members of each `Group N`, meaning the strains whose mismatches are exactly the same across all reads:

```
Group 0:
STRAIN_3
STRAIN_4
```

**`per_strain_llr.csv`** is written only with `-l`. It has the columns `names.arg,p,llr` for the top `-z` strains.

The working directory also keeps one `<STRAIN>.sam` file per reference strain from the most recent run.

## Memory footprint

Peak usage is roughly **12 bytes × reads × strains**. Grouping unidentifiable strains shrinks the matrices, so real usage is often lower. Use `-m` to keep runs within your machine's memory.

## Known limitations

- **`-m` is a hard limit.** If more aligned reads remain than `-m` allows, the program stops with `Error: Too many aligned reads!`.
- **Single-line FASTA only.** The MSA and the reference sequence files must have each sequence on one line. Wrapped sequences are read incorrectly.
- **SAM files must line up.** Reads are matched across references by line position, so every reference's SAM file must contain the same number of alignment lines. The program stops with an error if the counts differ.
- **Strain names can't start with `-`.** Such a name would be read as the next option.
- **`-n/--no-read-sam`** is accepted but does nothing yet.

## Repository layout

| File | Purpose |
|---|---|
| `eliminate_strains.c` | Entry point; runs the pipeline steps in order |
| `options.c` | Command-line parsing, including the multi-name `-r` option |
| `msa.c` | Reads the MSA |
| `align_references.c` | Locates reference strains in the MSA and loads their column maps |
| `bowtie_alignment.c` | Builds Bowtie2 indexes and runs alignments |
| `sam.c` | Reads SAM files and drops reads that aligned to no reference |
| `build_mismatch_matrix.c` | Counts each read's mismatches against every strain |
| `calculate_proportions.c` | Groups strains, builds likelihoods, runs EM, writes outputs |
| `external/covid_em.c` | EM update, log-likelihood and SQUAREM routines |
| `external/hashmap.c` | Hash map by David Leeds (MIT license), used to group strains |
| `file_utils.c` | File path helpers |
| `clean_reads.c` | Compiled but not currently called |

## License

GNU General Public License v3.0. See `LICENSE`.