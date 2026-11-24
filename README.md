# influenza_wastewater_surveillance

Estimates the relative proportions of influenza strains in a wastewater sequencing sample.

Reads are aligned to a few reference strains, projected onto a multiple sequence alignment (MSA) of known influenza genomes, and compared against every strain in the MSA. An expectation–maximization (EM) mixture model then estimates what proportion of the sample each strain makes up.

## How it works

`eliminate_strains` runs these steps in order:

1. **Read the MSA.** Loads every strain's aligned sequence and name.
2. **Place the reads on the MSA.** How this is done depends on the alignment mode (`-A`, see [Alignment modes](#alignment-modes)). The steps below describe the default mode, `references`.
3. **Locate the reference strains.** For each reference strain, reads `<reference-sequences-dir>/<STRAIN>.fasta`, finds that strain in the MSA, and loads its base-to-column map from the positions file (`-p`).
4. **Align reads with Bowtie2.** For each reference strain, builds a Bowtie2 index in `<bowtie2-indexes-dir>/<STRAIN>` if one doesn't exist yet, then aligns the reads with `bowtie2 --all`, writing `<working-dir>/<STRAIN>.sam`. Only the first (best) alignment of each read is used.
5. **Read the SAM files and drop unaligned reads.** A read (or read pair) is dropped only if it failed to align to every reference strain. For paired reads, a pair counts as aligned if either mate is mapped; when only one mate maps, that mate's bases are still used.
6. **Build the mismatch matrix.** For each read, its best alignment (fewest mismatches to its reference) is projected onto MSA columns, and its mismatches against every strain in the MSA are counted.
7. **Group indistinguishable strains.** Strains whose mismatch counts are identical across every read cannot be told apart by this sample, so they are merged into a single `Group N`.
8. **Estimate proportions with EM.** Each read's likelihood under strain *j* is

   L<sub>ij</sub> = ε<sup>m<sub>ij</sub></sup> · (1 − ε)<sup>b<sub>i</sub> − m<sub>ij</sub></sup>

   where m<sub>ij</sub> is the read's mismatch count against strain *j*, b<sub>i</sub> is its number of aligned bases, and ε is the error rate (`-e`). EM, accelerated with SQUAREM and parallelized with OpenMP, finds the strain proportions that maximize the total likelihood. It runs from two starting points, one informed by the data and one random, and prints how far apart the two answers are as a convergence check.
9. **Optional likelihood-ratio test** (`-l`). For the top `-z` strains, refits the model without each strain and reports a log-likelihood ratio measuring how much that strain is needed.

## Requirements

- **C compiler with OpenMP**, such as GCC.
- **zlib** (`-lz`). All input files may be gzip-compressed.
- **Bowtie2**. `bowtie2` and `bowtie2-build` must be on your `PATH`.
- **minimap2**, only for the alignment mode `minimap2`. `minimap2` must be on your `PATH`.

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
| Positions file | `-p` | Not needed in mode `anchor-imputed`. For each reference strain, a header line `>STRAIN` followed by one line of comma-separated integers: the MSA column index of each base of that strain's reference sequence, in order. The number of entries must equal the reference sequence length. |
| Reference sequences | `-g` | Not needed in mode `anchor-imputed`. A directory containing `<STRAIN>.fasta`, the non-imputed sequence, for every strain (members of redundant groups included). The header must be the strain name, and the sequence must be on **a single line**. |
| Bowtie2 indexes | `-b` | A directory for `<STRAIN>.1.bt2` and the other index files, and for the indexes of all strains that the modes `minimap2` (`all_reference_sequences.mmi`), `anchor-imputed` (`all_imputed_strains.*`) and `anchor-real` (`all_real_strains.*`) use. Missing indexes are built automatically, but the directory must already exist. |
| Redundant groups | `-G` | Optional. The groups file of a nonredundant MSA: a `>cN` line for each group, followed by the names of its member strains. Without it, strains that are members of a group cannot be chosen as reference strains. |
| Reads | `-1`/`-2` or `-0` | FASTQ by default, or FASTA with `-a`. Use `-P` with `-1` and `-2` for paired-end reads, or `-0` for single-end reads. |

## Usage

```bash
./eliminate_strains \
  -i msa.fasta \
  -p strain_positions.txt \
  -g reference_sequences \
  -b bowtie2_indexes \
  -r STRAIN_1 STRAIN_2 STRAIN_3 \
  -G redundant_groups.txt \
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
| `-p` | `--non-imputed-positions-filepath` | file | required (not in `anchor-imputed`) | Base-to-MSA-column map for each reference strain |
| `-g` | `--reference-sequences-dir` | dir | required (not in `anchor-imputed`) | Directory of `<STRAIN>.fasta`, the non-imputed sequence of every strain |
| `-b` | `--bowtie2-indexes-dir` | dir | required | Directory of Bowtie2 indexes, built if missing |
| `-A` | `--alignment-mode` | mode | `references` | How the reads are placed on the MSA: `references`, `minimap2`, `anchor-imputed` or `anchor-real` (see [Alignment modes](#alignment-modes)) |
| `-B` | `--bowtie2-args` | string | | Extra Bowtie2 options for aligning the reads, in every mode (quote them). They come after the program's own options |
| `-r` | `--reference-strains` | names | | Reference strains, space-separated. Required in mode `references`, optional in mode `minimap2`, not used in the anchor modes |
| `-G` | `--redundant-groups-filepath` | file | | Redundant groups file of the MSA, so that members of a group can be used as reference strains |
| `-I` | `--min-reference-identity` | decimal | 0.93 | Mode `minimap2`. A reference strain covers a read if at least 90% of the read aligns to it with at least this identity |
| `-M` | `--min-reads-per-added-reference` | int | 5 | Mode `minimap2`. A reference strain is added only if it covers at least this many reads that no reference strain covers yet |
| `-k` | `--min-strains-per-read` | int | 1 | Mode `minimap2`. A read is kept as an influenza read only if it aligns to at least this many strains, or to a strain given with `-r` |
| `-x` | `--minimap2-args` | string | see below | Mode `minimap2`. minimap2 options for aligning the reads to all strains, replacing the defaults (quote them; keep `-c`) |
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

## Alignment modes

`-A` chooses how the reads are placed on the MSA. Everything after that (mismatch matrix, grouping, EM) is the same in every mode.

| Mode | What aligns the reads | Reference strains |
|---|---|---|
| `references` (default) | Bowtie2, once per reference strain | `-r`, required |
| `minimap2` | minimap2 against all strains first, then Bowtie2 once per reference strain | `-r` optional; more are added from the reads |
| `anchor-imputed` | Bowtie2, once, against one index of every MSA row with its gaps removed | none |
| `anchor-real` | Bowtie2, once, against one index of the non-imputed sequence of every strain in `-g` | none |

- **`references`.** A read is used if it aligns to at least one reference strain. If it aligns to several, the reference it has the fewest mismatches to (per aligned base) places it.
- **`minimap2`.** See [Influenza reads and reference strains](#influenza-reads-and-reference-strains).
- **`anchor-imputed` and `anchor-real`.** No reference strains are chosen. Each read is aligned to all strains at once and only its best hit is kept (Bowtie2 without `--all`). That hit only says which MSA columns the read covers; the read is then compared with every MSA row on those columns, as in the other modes. In `anchor-imputed` a read can be placed by imputed bases, and base *i* of a row's sequence is the row's *i*-th column that is not a gap, so `-p` and `-g` are not needed. In `anchor-real` only real sequence can place a read, and the columns come from the positions file. `-r` is accepted and ignored in both, so that a script that always passes it can run every mode.
- **The index of all strains** is built on the first run (`all_imputed_strains.*` or `all_real_strains.*` in the `-b` directory) and reused. For the full MSA this takes a while and should be done by one run before any parallel runs. A small `.info` file next to it records what it was built from; a run with a different MSA (or a different number of strains in `-g`) stops with an error. Delete the index files to have it rebuilt.

`-B` adds Bowtie2 options in every mode, for example a looser score limit, so that more divergent reads align:

```bash
./eliminate_strains ... -r STRAIN_1 STRAIN_2 -B "--score-min L,-0.6,-1.0"
./eliminate_strains ... -A anchor-imputed
```

## Influenza reads and reference strains

This section describes the mode `minimap2` (`-A minimap2`). Before the Bowtie2 step, minimap2 aligns the reads to the non-imputed sequences of the strains (the `-g` directory). This does two things.

**It keeps only the influenza reads.** All reads are aligned to all strains. An alignment counts only if at least 90% of the read is aligned with an identity of at least 0.8. A read with such an alignment to at least `-k` strains, or to a strain given with `-r`, is an influenza read (in paired-end mode a pair is one read, and one mate is enough). The influenza reads are written to `influenza_reads_1.fastq` and `influenza_reads_2.fastq` (or `influenza_reads.fastq`) in the working directory, and the rest of the run, Bowtie2 included, uses only these.

**It completes the reference strains.** The strains given with `-r` are always used. A strain *covers* a read if at least 90% of the read aligns to it with an identity of at least `-I` (0.93). For the influenza reads that no given strain covers, strains are added:

1. All reads were aligned to all strains with limits that keep minimap2 fast: for each read it reports the best-scoring strain and up to 100 more that score close to it. From these lists, the program takes the strain that covers the most reads not covered yet, and repeats until every read is covered. This choice is provisional: a read that matches thousands of strains equally well was given only 101 of them, and two such reads need not have been given the same ones, so it has more strains than needed.
2. The reads are aligned again to the provisional strains only, with the limits lifted, so each read has the complete list of provisional strains that cover it.
3. The same greedy step on the complete lists gives the strains to add. It stops when the next strain would cover fewer than `-M` (5) more reads. Ties go to the strain with the most matching bases over all its alignments.

Without `-r`, all reference strains are chosen this way. The reference strains are listed in `selected_references.tsv` in the output directory, the given ones first.

Things to know:

- **An added strain never replaces a given one.** For each read, the program keeps the alignment to the reference strain with the lowest mismatch rate, so an added strain only changes the reads it fits better than the given ones.
- **The minimap2 index** of all sequences is built on the first run and kept as `all_reference_sequences.mmi` in the `-b` directory. Delete it after changing the `-g` directory. Let one run build it before starting runs in parallel.
- **Redundant groups.** A strain that is a member of a group has no MSA row of its own; its row is the group's (`cN`). With `-G`, such strains can be reference strains, and their mismatches are counted against the group's row at the strain's own positions. Without `-G`, they are ignored.
- **Strains without an MSA row are ignored**, also for deciding whether a read is an influenza read. To keep a problematic strain out, remove it from the MSA (and from the groups file).
- **minimap2 options.** All reads against all strains: `-x sr -c --secondary=yes -N 100 -p 0.9 -f 10000`, the options EsViritu aligns with; `-x "..."` replaces them (keep `-c`: without it the identity of an alignment is only approximate). Against the given strains and the provisional strains, the same options are followed by `-N 1000000 -f 1000000 -K 500k`.
- A read that is covered by no reference strain is still kept and still goes to Bowtie2, which may or may not align it.

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

**`selected_references.tsv`** lists the reference strains (the given ones, then the added ones in the order they were added), with each one's MSA row and the number of influenza reads it newly covers.

The working directory also keeps the influenza reads (`influenza_reads*.fastq`) and one `<STRAIN>.sam` file per reference strain from the most recent run.

## Memory footprint

Peak usage is roughly **12 bytes × reads × strains**. Grouping unidentifiable strains shrinks the matrices, so real usage is often lower. Use `-m` to keep runs within your machine's memory.

## Known limitations

- **`-m` is a hard limit.** If more aligned reads remain than `-m` allows, the program stops with `Error: Too many aligned reads!`.
- **Single-line FASTA only.** The MSA and the reference sequence files must have each sequence on one line. Wrapped sequences are read incorrectly.
- **Reads must have distinct names.** The lines of a read are recognised in a SAM file by its name; only the first line of each read (or of each mate) is used, so further alignments reported by `bowtie2 --all` are ignored.
- **Reads at the edge of a missing stretch** (modes `references`, `minimap2`, `anchor-real`). A non-imputed sequence simply skips what the strain is missing, so a read that starts inside a missing stretch can align across the join, and its first bases are then counted on the wrong MSA columns. `anchor-imputed` does not have this problem.
- **`anchor-imputed` can place a read by imputed bases.** If the imputation is wrong there, the read is still placed there.
- **Strain names can't start with `-`.** Such a name would be read as the next option.
- **`-n/--no-read-sam`** is accepted but does nothing yet.

## Repository layout

| File | Purpose |
|---|---|
| `eliminate_strains.c` | Entry point; runs the pipeline steps in order |
| `options.c` | Command-line parsing, including the multi-name `-r` option |
| `msa.c` | Reads the MSA |
| `align_references.c` | Locates reference strains in the MSA and loads their column maps |
| `select_references.c` | Mode `minimap2`: keeps the influenza reads and adds reference strains |
| `anchor_alignment.c` | Modes `anchor-imputed` and `anchor-real`: index of all strains, one alignment, MSA columns of each read's best hit |
| `redundant_groups.c` | Reads the redundant groups file |
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