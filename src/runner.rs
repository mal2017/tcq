use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use super::handler::Nascent;
use super::filter::*;
//use std::io::prelude::*;

pub fn run_through_bam(ib: &str, ob: &str, tag: &str, p: usize, blk: Option<&str>) {
	let filt: Option<ConvFilter> = match blk {
		Some(b) => {
			Some(ConvFilter::from_vcf_path(b, p).unwrap())
		},
		None => None,
	};
	info!("beginning run...");
	info!("opening bams...");
	// https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let hdr = bam::header::Header::from_template(bam.header());
	let mut obam = bam::Writer::from_path(ob, &hdr).unwrap();

	info!("setting thread usage...");

	if p >= 2 {
		let p2 = if (p % 2) == 0 {
			p / 2
		} else {
			(p - 1) / 2
		};
		bam.set_threads(p2).unwrap();
		obam.set_threads(p-p2).unwrap();
	} else {
		bam.set_threads(1).unwrap();
		obam.set_threads(1).unwrap();
	}

	info!("annotating reads with t>>c conversions...");
	// TODO make interval tree lookup for vcf filter
	// TODO optimize  mapping and pusing, maybe use mut iter
	/*bam.records().into_iter()
				.map(|a| a.unwrap())
				.map(|mut a| {a.push_tc_conv_aux(tag.as_bytes()).unwrap();a})
				.map(|a| obam.write(&a).unwrap())
				.for_each(drop);*/
}
