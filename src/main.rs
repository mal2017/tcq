#![forbid(unsafe_code)]

extern crate clap;
extern crate rust_htslib;
extern crate bio;
extern crate tcq;
#[macro_use] extern crate log;
extern crate env_logger;

use std::env::set_var;

fn main() {
	use clap::{Arg, App};
	use tcq::runner;
	use tcq::filter;



	let matches = App::new("tcq")
                          .version("0.1.1")
                          .author("Matt Lawlor <matt.a.lawlor@gmail.com>")
                          .about("T>>C conversion annotation util for working with nascent RNA species.\nAdds valid T>>C conversions to tag of your choice.\nRequires revcomp (-) seqs & MD tags.")
                          .arg(Arg::with_name("IBAM")
                               .help("bam to annotate")
                               .required(true)
                               .index(1))
                          .arg(Arg::with_name("OBAM")
                               .help("bam to write")
                               .required(true)
                               .index(2))
						  .arg(Arg::with_name("BLKLIST")
					  			.help("vcf, bcf, bed, or gff with blacklisted sites")
								.short("b")
								.long("blacklist")
								.takes_value(true))
                          .arg(Arg::with_name("TAG")
                          	   .help("tag to store tc conversion number")
                          	   .short("t")
                          	   .long("tag")
                          	   .takes_value(true))
						  .arg(Arg::with_name("CONTIG_PREFIX")
					  		   .help("prefix blacklist contigs with this for quick compatibility")
						   	   .short("c")
						   	   .long("contig-prefix")
						   	   .takes_value(true))
                          .arg(Arg::with_name("THREADS")
                          	   .help("threads to use")
                          	   .short("p")
                          	   .long("threads")
                          	   .takes_value(true))
                          .arg(Arg::with_name("VERBOSE")
                               .help("run in verbose mode")
                               .short("v")
                               .long("verbose"))
                          .get_matches();

    let verbose: bool = matches.is_present("VERBOSE");

    match verbose {
      true => {
        set_var("RUST_LOG","info");
      },
      false => {
        set_var("RUST_LOG","error")
      }
    }

    env_logger::init();

    info!("tcq successfully initialized...");

    let bam_file: &str = matches.value_of("IBAM").unwrap();
    let obam_file: &str = matches.value_of("OBAM").unwrap();
    let tag: &str = matches.value_of("TAG").unwrap_or("ZX");
    let threads: usize = matches.value_of("THREADS").unwrap_or("1").parse().unwrap();
	let blk: Option<&str> = matches.value_of("BLKLIST");
	let ctg_prefix: &str = matches.value_of("CONTIG_PREFIX").unwrap_or("");

    info!("arguments parsed...");
    // TODO: Check files are valid
	// TODO: mapq
	// TODO: per base phred cutoff
	// TODO: check md tag exists
    // TODO: Check tag starts with X, Y, or Z
	// TODO: Check tag not already in use
	// TODO: add force option if already in use
    // TODO: check revcomp

    runner::run_through_bam(bam_file, obam_file, tag, threads, blk, ctg_prefix);

    info!("tcq run complete");
}
