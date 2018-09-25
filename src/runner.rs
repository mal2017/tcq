use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use rayon::ThreadPoolBuilder;
use std::path::Path;
use std::str;
use bio::utils;
use super::handler::*;
use super::handler::Nascent;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
//use std::io::prelude::*;

pub fn run_through_bam(ib: &str, tag: &str, p: usize) {

	// TODO check reads have md tags,unless unmapped
	// https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let hdr = bam.header().to_owned();

	bam.set_threads(p);
	
	let mut bam_rec = bam::Record::new();
	let mut aux_md: String;
	let mut has_var: bool;
	let mut aux_nm: i64;

	//ThreadPoolBuilder::new().num_threads(2).build_global().unwrap();
	
	//let read_itr = bam.records().par_bridge();
	let read_itr = bam.records().into_iter();

	read_itr.map(|a| a.unwrap())
			.filter(|a| a.mutated())
			.map(|a| a.ref_seq())
			//.map(|a| {println!("{:?}",a);a})
			.for_each(drop);


}