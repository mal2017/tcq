use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use rayon::ThreadPoolBuilder;
use std::path::Path;
use std::str;
use bio::utils;
//use std::io::prelude::*;

pub fn run_through_bam(ib: &str, tag: &str, fa: &str, p: usize) {

	// TODO handle error
	let fa_path = Path::new(fa);

	// TODO handle error
	let mut fasta = fasta::IndexedReader::from_file(&fa_path);

	// TODO check reads have md tags,unless unmapped
	// https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let hdr = bam.header().to_owned();

	bam.set_threads(p);
	
	let mut bam_rec = bam::Record::new();
	let mut aux_md: String;
	let mut ref_seq: Vec<u8> = vec![b'a'];
	let mut tid: &[u8];

	while let Ok(_x) = bam.read(&mut bam_rec) {
		aux_md = String::from_utf8(bam_rec.aux(b"MD")
										  .unwrap()
										  .string()
										  .to_vec()).unwrap();

		//fasta.fetch("chr1",bam_rec.pos(),200).unwrap();
		//fasta.read(&mut ref_seq);

		
		//println!("{:?}", str::from_utf8(&ref_seq));
		/* Here, get helper that parses md tag
		things to consider:
			1. just get positions of mismatches that are 
			c's or g's (t>>c or a>>g)
			2. exclude those that are in deletions
		
		println!("{:?}", aux_md);
		*/
	}	
}