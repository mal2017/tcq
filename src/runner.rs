use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use super::handler::Nascent;
use super::handler::tid_2_contig;
use super::filter::*;
use super::spliced_read_utils::SplicedReadCigarStringView;
//use std::io::prelude::*;

pub fn run_through_bam(ib: &str, ob: &str, tag: &str, p: usize, blk: Option<&str>, pfx: &str, mq: u8) {
	let filt: Option<ConvFilter> = match blk {
		Some(b) => {
			Some(ConvFilter::from_vcf_path(b, p, pfx).unwrap())
		},
		None => None,
	};
	info!("beginning run...");
	info!("opening bams...");

	// https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let hdr = bam::header::Header::from_template(bam.header());
	let hdrv = bam.header().to_owned();
	let mut obam = bam::Writer::from_path(ob, &hdr).unwrap();

	let tid_lookup = tid_2_contig(&hdrv);

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
	bam.records().into_iter()
				.map(|a| a.unwrap())
				.filter(|a| a.mapq() >= mq )
				.map(|mut a| {a.push_tc_conv_aux(tag.as_bytes(), &filt, &tid_lookup).unwrap();a})
				.map(|a| obam.write(&a).unwrap())
				.for_each(drop);
}
