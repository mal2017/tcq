use bio::data_structures::interval_tree::IntervalTree;
use rust_htslib::bcf;
use rust_htslib::bam;
use std::collections::HashMap;
use bio::utils::Interval;
use rust_htslib::bcf::Read;

#[derive(Debug)]
pub struct ConvFilter {
    inner: Option<HashMap<String, IntervalTree<u32, u32>>>,
}

impl ConvFilter {
    pub fn from_vcf_path(v: &str, p: usize) -> Result<Self, ConvFilError> {
        let mut blank: HashMap<String, IntervalTree<u32, u32>> = HashMap::new();
        let mut vcf = bcf::Reader::from_path(v).unwrap();
        vcf.set_threads(p);

        Ok(
            ConvFilter {
                inner: None,
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
