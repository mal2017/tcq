use super::filter::*;
use super::handler::tid_2_contig;
use super::handler::LibraryType;
use super::handler::Nascent;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;

/// Iterates through bam records and annotates each read with the T>>C conversion count.
///
/// This this is the core of the tcq executable. It will write a bam to the  output
/// file specified by `ob`.
///
/// 1. ib: path for input bam with unannotated reads
/// 1. ob: output bam with annotated reads
/// 1. tag: preferred aux tag in which to store conversion count for each read
/// 1. p: threads to use for reading and writing
/// 1. blk: path to optional indexed vcf/bcf blacklist of individual positions to exclude
/// 1. mq: minimum read mapq
///
/// # Example (compiled, not run)
/// ```rust,no_run
/// use tcq::runner;
/// use tcq::handler::LibraryType::R1ANTISENSE;
/// runner::run_through_bam("in.bam", "out.bam", "XZ", 4, Some("filt.bcf"), 30, R1ANTISENSE);
/// ```
pub fn run_through_bam(
    ib: &str,
    ob: &str,
    tag: &str,
    p: usize,
    blk: Option<&str>,
    mq: u8,
    library: LibraryType,
) {
    info!("beginning run...");

    // Creates either a working filter or a none.
    let filt: Option<ConvFilter> = match blk {
        Some(b) => {
            info!("creating filter...");
            Some(ConvFilter::from_vcf_path(b, p).unwrap())
        }
        None => None,
    };

    info!("opening bams...");

    // Make bam reader
    let mut bam = bam::Reader::from_path(ib).unwrap();

    // Make bam header for writing output and seqname/target id conversion.
    let hdr = bam::header::Header::from_template(bam.header());
    let hdrv = bam.header().to_owned();

    // Initialize bam writer.
    let mut obam = bam::Writer::from_path(ob, &hdr).unwrap();

    // Create lookup hash table for converting TID to human readable chrom name.
    let tid_lookup = tid_2_contig(&hdrv);

    info!("setting thread usage...");

    // Set thread usage for reading/writing bam.
    if p >= 2 {
        let p2 = if (p % 2) == 0 { p / 2 } else { (p - 1) / 2 };
        bam.set_threads(p2).unwrap();
        obam.set_threads(p - p2).unwrap();
    } else {
        bam.set_threads(1).unwrap();
        obam.set_threads(1).unwrap();
    }

    // Begin iterator chain that processes each read individually.
    // 1-unwrap read
    // 2-make sure it is mapped
    // 3-confirm mapq > cutoff
    // 4-annotate with conversion count (see tcq::handler)
    // 5-write read to bam
    info!("annotating reads with t>>c conversions...");
    bam.records()
        .into_iter()
        .map(|a| a.unwrap())
        .filter(|a| !a.is_unmapped())
        .filter(|a| a.mapq() >= mq)
        .map(|mut a| {
            a.push_tc_conv_aux(tag.as_bytes(), &filt, &tid_lookup, &library)
                .unwrap();
            a
        }).map(|a| obam.write(&a).unwrap())
        .for_each(drop);
}
