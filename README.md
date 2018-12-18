# tcq

`tcq` is a utility written in the rust language for preprocessing of aligned NGS reads generated from
s<sup>4</sup>U labelled RNAseq experiments. `tcq` will count T to C base
conversions introduced as part of the workflow and add this to the auxiliary sam/bam tag
of your choosing. For more details on the protocol
itself see [Timelapse-seq](https://doi.org/10.1038/nmeth.4582) and
[SLAM-seq](https://doi.org/10.1126/science.aao2793).

## Motivation

Labelled nascent RNA-seq assays are perfect for measuring highly kinetic cellular
responses to chemical perturbation. To my knowledge, no purpose-built reusable tools
exist as of yet for preprocessing of this sort of data. It's not too hard to write
python or R scripts that do similar things, but I wanted to learn rust because it is
extremely fast, sort of fun, and can compile to many targets, including 100% static MUSL.

## Features

* Annotate individual reads with the number of s<sup>4</sup>U induced base conversions.
* Library-specific processing: first strand, second strand, and unstranded.
* Filter individual locations with a bcf/vcf file.
* Filter reads by mapq.
* Simple CLI.
* Multithreaded I/O.

## Getting Started

```bash
git clone https://github.com/mal2017/tcq.git

cd tcq

cargo install
```

On macOS mojave you may need to:
```bash
CC=clang cargo install
```

```bash
tcq -p 4 --tag ZX -m 30 --blacklist all_TC-AG_SNPs.bcf --r1-sense in.bam out.bam
```

I have plans to release static binaries, a conda package, and a docker image as the
project stabilizes.

## Running tests

```bash
cd tcq/

cargo test
```

## Credits

* https://github.com/rust-bio
* https://github.com/samtools/htslib
