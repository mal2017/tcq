

#[macro_use]
extern crate criterion;

use criterion::Criterion;
use criterion::black_box;

extern crate tcq;
extern crate rust_htslib;


use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::path::Path;
use tcq::handler::*;


//use handler::{LibraryType, Nascent};


fn forward_read() {
    let lib = LibraryType::R1SENSE;
    let bampath = Path::new("test/insertion_forward.r1s.myc.bam");
    let mut bam = bam::Reader::from_path(bampath).unwrap();
    let hdrv = bam.header().to_owned();
    let tid_lookup = tid_2_contig(&hdrv);
    let tcc: Vec<u32> = bam
        .records()
        .map(|a| a.unwrap())
        .into_iter()
        .map(|a| a.tc_conversions(&None, &tid_lookup, &lib))
        .collect();

}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("forward_read()", |b| b.iter(|| forward_read()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
