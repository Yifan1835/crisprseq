# CRISPRDecode paired-guide test data

**New users:** start with the [step-by-step tutorial](../../docs/usage/crisprdecode_tutorial.md). It includes a small first run and a 6-million-read-pair, 10,000-construct benchmark with exact truth checks.

These files are a fully synthetic, distributable truth set created for the nf-core/crisprseq test suite. They do not contain collaborator or biological sequencing data.
They can be reproduced with `python3 tests/crisprdecode/generate_fixtures.py`.

The paired FASTQ fixture contains five synchronized read pairs:

- two uniquely assignable signatures (`construct_a` and `construct_b`);
- one signature shared by `construct_dup_1` and `construct_dup_2`;
- one extracted signature absent from the library;
- one R1 sequence shorter than the required five-base extraction window.

`construct_zero` is uniquely identifiable but deliberately unobserved. The expected assignment totals are therefore two unique, one ambiguous, one unassigned and one extraction failure.

Run the Python tests with:

```bash
python3 -m unittest discover -s tests/crisprdecode -v
```

For the fixed collaborator test version, GitHub Actions instructions and feedback template, see [COLLABORATOR_TESTING.md](COLLABORATOR_TESTING.md).
