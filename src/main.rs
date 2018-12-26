#![forbid(unsafe_code)]

extern crate bio;
extern crate clap;
extern crate rust_htslib;
extern crate tcq;
#[macro_use]
extern crate log;
extern crate env_logger;

use std::env::set_var;
use tcq::validators::*;

fn main() {
    use clap::{App, Arg, ArgGroup};
    use tcq::handler;
    use tcq::runner;

    let matches = App::new("tcq")
                          .version("0.2.0")
                          .author("Matt Lawlor <matt.a.lawlor@gmail.com>")
                          .about("Util for SLAM-/Timelapse-seq. Adds valid T>>C conversions to tag of your choice.\nRequires revcomp (-) seqs & MD tags.")
                          .arg(Arg::with_name("IBAM")
                               .help("unannotated reads")
                               .required(true)
                               .index(1)
						   	   .validator(bam_seems_ok))
                          .arg(Arg::with_name("OBAM")
                               .help("path to write annotated bam reads")
                               .required(true)
                               .index(2)
							   .validator(dir_exists))
					      .arg(Arg::with_name("R1SENSE")
					  		    .long("r1-sense"))
						  .arg(Arg::with_name("R1ANTISENSE")
					  			.long("r1-antisense"))
						  .arg(Arg::with_name("UNSTRANDED")
					  			.long("unstranded"))
						  .group(ArgGroup::with_name("LIBRARY")
					  			.args(&["R1SENSE","R1ANTISENSE","UNSTRANDED"])
								.required(true))
						  .arg(Arg::with_name("BLKLIST")
					  			.help("indexed vcf or bcf with blacklisted sites")
								.short("b")
								.long("blacklist")
								.takes_value(true)
								.validator(vcf_seems_ok))
                          .arg(Arg::with_name("TAG")
                          	   .help("tag to store tc conversion number")
                          	   .short("t")
                          	   .long("tag")
                          	   .takes_value(true)
							   .validator(tag_is_reserved_local))
                          .arg(Arg::with_name("SOFTCLIPS")
                               .help("use if softclips are present in your alignments")
                               .short("s")
                               .long("softclips"))
						  .arg(Arg::with_name("MAPQ")
					  		   .help("MAPQ must be greater than OR EQUAL TO provided cutoff; default 0")
						   	   .long("mapq")
						   	   .short("m")
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
            set_var("RUST_LOG", "info");
        }
        false => set_var("RUST_LOG", "error"),
    }

    env_logger::init();

    info!("tcq successfully initialized...");

    let bam_file: &str = matches.value_of("IBAM").unwrap();
    let obam_file: &str = matches.value_of("OBAM").unwrap();
    let tag: &str = matches.value_of("TAG").unwrap_or("ZX");
    let threads: usize = matches.value_of("THREADS").unwrap_or("1").parse().unwrap();
    let mapq: u8 = matches.value_of("MAPQ").unwrap_or("0").parse().unwrap();
    let blk: Option<&str> = matches.value_of("BLKLIST");
    let softclips: bool = matches.is_present("SOFTCLIPS");

    let library = if matches.is_present("R1SENSE") {
        handler::LibraryType::R1SENSE
    } else if matches.is_present("R1ANTISENSE") {
        handler::LibraryType::R1ANTISENSE
    } else {
        handler::LibraryType::UNSTRANDED
    };

    info!("arguments parsed...");

    // TODO: per base phred cutoff
    // TODO: check md tag exists
    // TODO: add force option if already in use
    // TODO: check revcomp

    runner::run_through_bam(bam_file, obam_file, tag, threads, blk, mapq, library, softclips);

    info!("tcq run complete");
}
