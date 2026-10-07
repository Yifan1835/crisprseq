# CRISPRDecode: first run and a larger test

This guide starts with synthetic data, so you do not need to prepare your own
sequencing files. Run commands in a terminal from the root of the **CRISPRDecode
test checkout** (the folder containing `main.nf` and `bin/crisprdecode_assign.py`).

CRISPRDecode takes two synchronized reads from each DNA fragment, extracts one
guide from R1 and one from R2, and counts the matching **construct**. A construct
is a library row describing a guide pair. Both guides must match exactly.
This version supports paired-end, paired-guide libraries without UMI/iBAR.
It does not perform UMI deduplication or prime-editing library decoding.

## 1. Start from the supplied code

Ask the project owner for the checkout/archive that contains this guide and
`tests/crisprdecode/benchmark_scale.py`. The new tutorial and scale benchmark
must be included in the version you receive; older test releases may not have them.

```bash
cd /path/to/crisprseq
git rev-parse HEAD
git status --short
python3 --version
ls bin/crisprdecode_assign.py tests/crisprdecode/benchmark_scale.py
```

Replace `/path/to/crisprseq` with your actual folder. Record the commit SHA and
any local changes. An archive without Git metadata should be unpacked into a
Git checkout before running this benchmark, which records the commit automatically.
Do not assume that `nextflow run nf-core/crisprseq` downloads this test feature.
All Nextflow examples below run the code in your current directory using `.`.

There are two execution levels:

| Level | Requirements | What it checks |
| --- | --- | --- |
| Python decoding | Python 3.12.8 used for this test; standard library only | Library validation, real FASTQ extraction/assignment, multi-sample aggregation, exact truth comparison |
| Nextflow pipeline | Java 17, Nextflow 25.04.0, Docker; nf-test 0.9.3 for developer tests | Workflow execution and process/container integration |

The Nextflow versions above match the existing
[collaborator CI configuration](../../.github/workflows/crisprdecode-collaborator.yml).
Python decoding is the easiest first run. It does not require Conda packages,
Docker, Nextflow, a reference genome, or network access.

## 2. Run a small example first

```bash
python3 -m unittest discover -s tests/crisprdecode -v
python3 tests/crisprdecode/benchmark_scale.py \
  --outdir /tmp/crisprdecode-first-run \
  --constructs 100 --samples 2 --pairs-per-sample 1000
```

Choose a new output directory if that path already exists. The benchmark refuses
to overwrite a previous run. The existing unit suite should pass 8 tests, and the
benchmark should finish with `"status": "PASS"`.

The benchmark generates paired compressed FASTQs, a library, a samplesheet, and
expected results. It then runs the three actual Python decoding scripts in order:
validate the library → assign each sample → combine the counts. It compares the
complete count matrix, assignment summary, and library recovery table against
truth generated from the simulated read identities, independently of the decoder.

Open these files:

```bash
cat /tmp/crisprdecode-first-run/results/assignment_summary.tsv
head -n 6 /tmp/crisprdecode-first-run/results/count_table.count.txt
cat /tmp/crisprdecode-first-run/report.json
```

For **each** of the two samples you should see:

| Metric | Expected read pairs |
| --- | ---: |
| total_reads | 1,000 |
| extracted_reads | 950 |
| unique_reads | 850 |
| ambiguous_reads | 50 |
| unassigned_reads | 50 |
| extraction_failed_reads | 50 |

Here, `total_reads` counts **pairs**, not individual R1 plus R2 records.
`extracted_reads` overlaps the assignment categories; do not add it to the total.
The checks are `850 + 50 + 50 + 50 = 1,000` and `850 + 50 + 50 = 950`.

## 3. Run the medium/large synthetic test

Use a disk location with space for generated FASTQs and results. The default is
10,000 library rows, 6 samples, and 1,000,000 read pairs per sample: **6 million
pairs / 12 million FASTQ records** in total. Samples run sequentially in the
Python benchmark; it measures local decoding rather than cluster throughput.

```bash
python3 tests/crisprdecode/benchmark_scale.py \
  --outdir /tmp/crisprdecode-scale \
  --constructs 10000 --samples 6 --pairs-per-sample 1000000
```

Expected per sample: 850,000 unique, 50,000 ambiguous, 50,000 unassigned, and
50,000 extraction failures. The count matrix has 10,000 data rows and 8 columns:
`sgRNA`, `Gene`, then 6 sample columns. Despite the inherited column name
`sgRNA`, each row represents a **construct**, and `Gene` contains its `target_id`.

The library contains 8,998 observable unique constructs, two rows sharing one
signature, and 1,000 deliberately unobserved unique constructs. Shared individual
guide sequences would be allowed; ambiguity is determined by the **pair**.
Duplicate-signature rows receive zero counts. Random sampling can leave additional
unique rows unobserved at smaller read depths; the generated truth accounts for this.

The simulated reads are normally 150 bases long. Failure cases include missing
anchors and short extraction windows. R2 is reverse-complemented in the FASTQ.
The category mixture is fixed, and seed 20261007 makes guide selection reproducible.
Control/treatment labels are synthetic sample labels; no biological effect is
simulated and no hit-recovery claim follows from this test.

Output layout:

```text
crisprdecode-scale/
├── inputs/       # library, samplesheet, FASTQ.gz files, Nextflow params.json
├── truth/        # independently generated expected tables
├── results/      # actual Python outputs
├── logs/         # stdout/stderr of each decoding command
└── report.json   # checks, commands, timings, Python/OS, commit and SHA-256 hashes
```

The report records generation time separately from validation, each sample's
decoding time, and aggregation time. Total time includes hashing and verification.
It does not measure peak memory.

## 4. Optional: run the same data through Nextflow

Complete step 2 or 3 first. The generated `inputs/params.json` contains absolute
input paths, the correct extraction parameters, and a separate Nextflow output
directory. If you copy the dataset to another computer, update these paths and
the paths in `samplesheet.csv`.

Check the environment:

```bash
java -version
nextflow -version
docker info
```

Run from the same test checkout as step 1:

```bash
NXF_VER=25.04.0 nextflow run . \
  -profile docker \
  -params-file /tmp/crisprdecode-scale/inputs/params.json \
  -work-dir /tmp/crisprdecode-scale/nextflow_work
```

Do not use `-stub` for this data test: stub processes do not decode reads.
Do not set `--bowtie`, trimming adapters, UMI options, or downstream statistical
flags for this first comparison. Extraction parameters describe the reads reaching
the decoder; extra trimming would change their geometry.

After Nextflow completes, verify its published outputs with the same truth:

```bash
python3 tests/crisprdecode/benchmark_scale.py \
  --outdir /tmp/crisprdecode-scale \
  --verify-results /tmp/crisprdecode-scale/nextflow_results/crisprdecode
```

All four checks must report `PASS`. See
[output documentation](../output/screening.md#crisprdecode-paired-guide-count)
for the published files. To resume an interrupted Nextflow run, repeat its command
with `-resume`, retaining the same work directory and launch directory.

This is a runnable procedure to validate Nextflow where its dependencies are
available. A Python benchmark PASS alone does **not** verify Nextflow, Docker,
MultiQC, MAGeCK statistics, or real sequencing data. The separate
[developer/CI guide](../../tests/crisprdecode/COLLABORATOR_TESTING.md) explains
nf-test coverage; its screening-routing tests use stubs.

## 5. Replace the synthetic inputs with your own data

Prepare two files, following the examples in the generated `inputs/` directory.

**Samplesheet: CSV, with this header and one row per sample.**

```csv
sample,fastq_1,fastq_2,condition
control_1,/data/control_1_R1.fastq.gz,/data/control_1_R2.fastq.gz,control
treated_1,/data/treated_1_R1.fastq.gz,/data/treated_1_R2.fastq.gz,treatment
```

Use unique sample names without spaces and synchronized R1/R2 files. Both read
files are required for CRISPRDecode. Use absolute paths.

**Construct library: TSV, with a header and actual tab separators.**

```tsv
construct_id	target_id	spacer_r1	spacer_r2
construct_A	TARGET_A	AAAAA	CCCCC
construct_B	TARGET_B	GGGGG	TTTTT
```

The five-base guides here illustrate the file format; use your actual guide
sequences. IDs must be unique and all fields non-empty. Sequences accept A/C/G/T
(case is normalized to uppercase). R1 spacers must all have the same length;
R2 spacers must all have the same length, which can differ from R1's length.
Keep one target label per construct. This is not the three-column MAGeCK library
format described for the default route.

**Set extraction parameters from your actual read structure.**

The decoder optionally reverse-complements the whole read **first**, then finds
the first exact anchor, skips the offset bases after it, and extracts the number
of bases specified by that read end's library spacer length.

For the illustrative library above:

```text
R1 as seen by decoder: ACGTAC GG AAAAA TTT...
                      anchor 2  guide tail
R2 after reverse complement: TGCATG GG CCCCC TTT...
                             anchor 2  guide tail
```

This requires R1 anchor `ACGTAC`, R2 anchor `TGCATG`, offsets 2, and R2 reverse
complement enabled, just as in the synthetic benchmark. If a read starts directly
with the guide, omit its anchor and use offset 0. Without an anchor, offset is
the zero-based start position in the oriented read. Do not copy example anchors
into a real run without checking the sequencing design.

A real-data command, for this particular geometry, is:

```bash
NXF_VER=25.04.0 nextflow run . -profile docker \
  --analysis screening --screening_count_method crisprdecode \
  --input /data/samplesheet.csv --library /data/construct_library.tsv \
  --crisprdecode_r1_anchor ACGTAC --crisprdecode_r2_anchor TGCATG \
  --crisprdecode_r1_offset 2 --crisprdecode_r2_offset 2 \
  --crisprdecode_reverse_complement_r2 \
  --outdir /data/crisprdecode_results
```

First inspect assignment QC and counts. Configure contrasts and downstream
statistics separately using the [screening guide](screening.md).
An output count table alone is not a biological hit list.

## 6. Understand unexpected results

| Symptom | Meaning / first check |
| --- | --- |
| Missing benchmark script or unknown CRISPRDecode parameter | Confirm the supplied checkout contains the feature; record its commit and local changes. |
| Output directory already exists | Use a new directory; do not delete an earlier result to rerun. |
| Duplicate column / missing column / nonuniform spacer length | Check the TSV header, tab separators, and library lengths. |
| High extraction_failed_reads | Check anchors, offsets, orientation, and read length. Either end can fail. |
| High unassigned_reads | Extraction succeeded, but the pair is absent from the library; check library version, guide orientation, and sequencing mismatches. |
| ambiguous_reads | The extracted pair matches multiple library rows; no row receives those counts. |
| Zero-count construct | Inspect assignability: the construct may be unobserved or have an indistinguishable signature. |
| Desynchronised reads / different numbers of records | Pair names or record counts disagree. Repair input pairing; do not bypass the check. |
| Docker or plugin download error | Save the first error and environment versions; decoding accuracy has not yet been tested. |
| Truth mismatch | Preserve actual and expected files and logs. Do not edit truth or accept snapshots to force a pass. |

## 7. Send useful feedback

Copy this short form to your reply to the project owner:

```text
Commit SHA and local changes:
Step number / first sentence that was unclear:
Command I ran:
Expected result:
Actual result / first error:
Python benchmark: PASS / FAIL / NOT RUN
Nextflow plus truth check: PASS / FAIL / NOT RUN
Path to report.json and logs (or CI run URL):
Setup time:
Suggested wording:
```

Feedback about a confusing instruction is useful even if no command has failed.
