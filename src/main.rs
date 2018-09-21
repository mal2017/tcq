#![forbid(unsafe_code)]

#[macro_use]
extern crate clap;
extern crate rust_htslib;
extern crate bio;
extern crate tcq;


fn main() {
	use clap::{Arg, App, SubCommand};
	use tcq::runner;

	let matches = App::new("tcq")
                          .version("0.0.9999")
                          .author("Matt Lawlor <matt.a.lawlor@gmail.com>")
                          .about("T>>C conversion util")
                          .arg(Arg::with_name("BAM")
                               .help("bam to use")
                               .required(true)
                               .index(1))
                          .arg(Arg::with_name("TAG")
                          	   .help("tag to store tc conversion number")
                          	   .short("t")
                          	   .long("tag")
                          	   .takes_value(true))
                          .arg(Arg::with_name("THREADS")
                          	   .help("threads to use")
                          	   .short("p")
                          	   .long("threads")
                          	   .takes_value(true))
                          .get_matches();

    let bam_file: &str = matches.value_of("BAM").unwrap();
    let tag: &str = matches.value_of("TAG").unwrap_or("ZX");
    let threads: usize = matches.value_of("THREADS").unwrap_or("1").parse().unwrap();

    // TODO: Check files are valid
    // TODO: Check fasta is indexed
    // TODO: Check tag starts with X, Y, or Z

    runner::run_through_bam(bam_file, tag, threads);
    
}
