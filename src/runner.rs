use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use std::path::Path;
use std::str;
use bio::utils;
use super::handler::*;
use super::handler::Nascent;
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

	let read_itr = bam.records().into_iter();

	// TODO check bits for revcomp set on - strand reads
	// TODO make interval tree lookup for vcf filter
	read_itr.map(|a| a.unwrap())
			.filter(|a| a.mutated())
			.map(|a| a.get_tc_mismatch_pos())
			//.map(|a| {println!("{:?}",a);a})
			.for_each(drop);


}