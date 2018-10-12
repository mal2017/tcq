use bio::data_structures::interval_tree::IntervalTree;
use rust_htslib::bcf;
use rust_htslib::bam;
use std::collections::HashMap;
use rust_htslib::bcf::Read;
use std::str;
use core::ops::Range;

#[derive(Debug)]
pub struct ConvFilter {
    pub inner: HashMap<String, IntervalTree<u32, u32>>,
}

impl ConvFilter {
    pub fn from_vcf_path(v: &str, p: usize, pfx: &str) -> Result<Self, ConvFilError> {
        info!("creating blacklist filter...");
        let mut blank: HashMap<String, IntervalTree<u32, u32>> = HashMap::with_capacity(100);
        let mut vcf = bcf::Reader::from_path(v).unwrap();
        vcf.set_threads(p).unwrap();
        let hdr = vcf.header().to_owned();

        let mut vcf_records = vcf.records();

        let mut record: bcf::record::Record;
        let mut pos: u32;
        //let mut rng: Range<u32>;
        let mut chrom: String;

        // https://doc.rust-lang.org/std/collections/hash_map/enum.Entry.html
        while let Some(r) =  vcf_records.next() {
            record = r.unwrap();
            chrom = format!("{}{}",pfx,str::from_utf8(hdr.rid2name(record.rid().unwrap())).unwrap().to_owned());
            //println!("{:?}", chrom);
            pos = record.pos();
            blank.entry(chrom)
                     .and_modify(|a| a.insert(Range  {start: pos, end: pos+1},0))
                     .or_insert({ let mut a = IntervalTree::new();
                                  a.insert(Range  {start: pos, end: pos+1},0);
                                  a
                                  });
        };
        Ok(
            ConvFilter {
                inner: blank,
            }
        )
    }
}


quick_error! {
    #[derive(Debug, Clone)]
    pub enum ConvFilError {
        NewFiltFromVcfError {
            description("Failed to make filter")
        }
    }
}

pub trait ToIntTree {
    fn add(&self, &mut ConvFilter) -> Result<(),ConvFilError> {
        Ok(())
    }
}

pub trait Blacklisted {
    fn is_blacklisted(&self) -> bool;
}



impl Blacklisted for bam::Record {
    fn is_blacklisted(&self) -> bool {
        false
    }
}

// generalized reader

// record to
