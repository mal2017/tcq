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

	bam.set_threads(1);
	
	let mut bam_rec = bam::Record::new();
	let mut aux_md: String;
	let mut has_var: bool;
	let mut aux_nm: i64;

	ThreadPoolBuilder::new().num_threads(p-1).build_global().unwrap();
	
	let read_itr = bam.records().par_bridge();

	read_itr.map(|a| a.unwrap())
			.filter(|a| a.mutated())
			.map(|a| a.ref_seq())
			.for_each(drop);

	
	//while let Ok(_x) = bam.read(&mut bam_rec) {
		//aux_md = String::from_utf8(bam_rec.aux(b"MD")
		//								  .unwrap()
		//								  .string()
		//								  .to_owned())
		//								  .unwrap();

		//aux_nm = bam_rec.aux(b"NM").unwrap().integer();

		//println!("{:?}", aux_nm);
		//has_var = handler::has_var(&aux_md);
		//println!("{:?} : {:?}", aux_md, has_var);

		//println!("{:?}", aux_md);
		
	//}	
}