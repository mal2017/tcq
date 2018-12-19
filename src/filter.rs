use bio::data_structures::interval_tree::IntervalTree;
use core::ops::Range;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::collections::HashMap;
use std::str;

/// This class is a collection of IntervalTrees wrapped into a HashMap.
#[derive(Debug)]
pub struct ConvFilter {
    /// Each entry in the HashMap represents a contig.
    pub inner: HashMap<String, IntervalTree<u32, u32>>,
}

impl ConvFilter {
    /// Create a new ConvFilter from a vcf/bcf file path.
    pub fn from_vcf_path(v: &str, p: usize) -> Result<Self, ConvFilError> {
        // Initialize empty hashmap holding interval trees.
        let mut blank: HashMap<String, IntervalTree<u32, u32>> = HashMap::with_capacity(100);

        // Open vcf.
        let mut vcf = bcf::Reader::from_path(v).unwrap();
        vcf.set_threads(p).unwrap();

        // Grab header for each TID/contig mapping.
        let hdr = vcf.header().to_owned();

        // Create iterator from records in vcf/bcf.
        let mut vcf_records = vcf.records();

        // Initialize variable to hold inidividual record so no reallocation necessary.
        let mut record: bcf::record::Record;
        let mut pos: u32;
        let mut chrom: String;

        // Iterate over records, storing them in `record` for each iteration of while loop.
        while let Some(r) = vcf_records.next() {
            record = r.unwrap();

            // Get chromosome name.
            chrom = format!(
                "{}",
                str::from_utf8(hdr.rid2name(record.rid().unwrap()))
                    .unwrap()
                    .to_owned()
            );

            // Get position of record.
            pos = record.pos();

            // Add an entry to the HashMap<IntervalTree> built above.
            blank
                .entry(chrom)
                .and_modify(|a| {
                    a.insert(
                        Range {
                            start: pos,
                            end: pos + 1,
                        },
                        0,
                    )
                }).or_insert({
                    let mut a = IntervalTree::new();
                    a.insert(
                        Range {
                            start: pos,
                            end: pos + 1,
                        },
                        0,
                    );
                    a
                });
        }
        Ok(ConvFilter { inner: blank })
    }
}

quick_error! {
    /// Errors if filter creation fails due to i/o or other file issue.
    #[derive(Debug, Clone)]
    pub enum ConvFilError {
        NewFiltFromVcfError {
            description("Failed to make filter")
        }
    }
}
