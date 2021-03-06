use super::expander::md_expanded;
use super::filter::ConvFilter;
use super::spliced_read_utils::*;
use regex::Regex;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::HeaderView;
use std::collections::HashMap;
use std::ops::Range;
use std::str;

/// Main engine for T>>C annotation.
pub trait Nascent {
    fn is_possible_nascent(&self) -> bool;

    fn md_ref_seq(&self) -> String;

    fn md_tag(&self) -> String;

    fn cand_tc_mismatch_ref_pos(&self) -> Vec<i32>;

    fn cand_tc_mismatch_pos_tuples(&self) -> Vec<(i32, i32)>;

    fn tc_conversion_pos(
        &self,
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>
    ) -> Vec<((i32, i32), bool)>;

    fn tc_conversions(
        &self,
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>,
    ) -> u32;

    fn push_tc_conv_aux(
        &mut self,
        auxtag: &[u8],
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>
    ) -> Result<(), NascentMolError>;
}

/// Trait for bam::Record that implemements methods for checking T>C and A>G conversion rates.
impl Nascent for Record {
    /// Initial filter that checks MD tag to see if an A or T mismatch exists.
    fn is_possible_nascent(&self) -> bool {
        // TODO check reads have md tags,unless unmapped
        if self.is_unmapped() {
            return false;
        };

        lazy_static! { // for speeeeed
            static ref foo_re: Regex  = Regex::new(r"A").unwrap();
            static ref bar_re: Regex  = Regex::new(r"T").unwrap();
        }

        let md = self.md_tag();

        foo_re.is_match(&md) | bar_re.is_match(&md)
    }

    /// Pull MD tag from record.
    fn md_tag(&self) -> String {
        String::from_utf8(self.aux(b"MD").unwrap().string().to_owned()).unwrap()
    }

    /// Expands the MD tag to a human readable string.
    ///
    /// See tcq::expander for details.
    fn md_ref_seq(&self) -> String {
        md_expanded(self.md_tag())
    }

    /// Finds candidate T>C positions relative to reference.
    fn cand_tc_mismatch_ref_pos(&self) -> Vec<i32> {
        // Get expanded MD tag
        let exp_md = self.md_ref_seq();

        // Create regex objects for what the reference base will be for a given
        // converted base. This is specific to the library orientation.
        lazy_static! { // for speeeeed
            static ref forward_re: Regex  = Regex::new(r"T").unwrap();
            static ref reverse_re: Regex  = Regex::new(r"A").unwrap();
        }

        match self.is_reverse() {
            true => reverse_re.find_iter(&exp_md)
                           .map(|a| a.start() as i32)
                           .collect(),
            false => forward_re.find_iter(&exp_md)
                           .map(|a| a.start() as i32)
                           .collect(),
        }

    }

    /// Makes a vector of tuples containing read and reference candidate T>>C positions.
    fn cand_tc_mismatch_pos_tuples(&self) -> Vec<(i32, i32)> {
        let cig = self.cigar();
        let start = self.pos() as i32;

        // from md tag so no softclips incl in ref pos,
        let cand_ref_pos = self.cand_tc_mismatch_ref_pos();
        // Make tuple with form (read, ref-pos)
        cand_ref_pos
            .iter()
            // use query pos calculated fron POS + (position within md tag)
            // use_softclips/use_dels is often irrelevant, all positions generally come From
            // highly covered variant calls, but good to set softclips to true just in case,
            // and include_dels to false so we don't get spurious tc calls.
            .map(|a| cig.read_pos_spliced(start + a, true, false, start))
            .map(|a| a.unwrap().unwrap_or(-1))
            .zip(cand_ref_pos.clone().into_iter())
            .filter(|a| a.0 >= 0)
            .collect()
    }

    /// Adds a bool to vector of read/ref position tuples.
    fn tc_conversion_pos(
        &self,
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>
    ) -> Vec<((i32, i32), bool)> {
        let mut cand_pos_tuples = self.cand_tc_mismatch_pos_tuples();
        let mut rng: Range<u32> = Range { start: 0, end: 0 };
        let chr = h.get(&(self.tid() as u32)).unwrap();
        let refpos = self.pos() as i32;
        let it;

        // BELOW THIS IS THE FILTERING PORTION
        cand_pos_tuples = match f {
            Some(k) => {
                it = &k.inner;
                cand_pos_tuples
                    .into_iter()
                    .filter(|a| !{
                        // negate false if overlap w. filt returns anything
                        rng = Range {
                            start: (refpos + a.1) as u32, // equiv to POS + (MD derived pos within read(Ie no softclips))
                            end: (refpos + a.1 + 1) as u32,
                        };
                        match it.get(chr) {
                            Some(j) => j.find(&rng).any(|_a| true),
                            None => false,
                        }
                    })
                    .collect()
            }
            None => cand_pos_tuples,
        };

        //println!("{:?}", cand_pos_tuples);
        // ABOVE THIS IS THE FILTERING PORTION
        let read_seq = self.seq();

        let is_r1 = self.is_first_in_template() || !(self.is_paired());

        let conv_target: Vec<u8> = match self.is_reverse() {
                true => vec![b'G'],
                false => vec![b'C'],
        };

        let enc_base_hit_itr = cand_pos_tuples
            .iter()
            .map(|a| a.0 as usize) // get read pos
            .map(|a| read_seq.as_bytes()[a])
            .map(|a| conv_target.contains(&a));

        cand_pos_tuples
            .clone()
            .into_iter()
            .zip(enc_base_hit_itr)
            .filter(|a| a.1)
            .collect()
    }

    /// Counts T>>C positions passing all filters.
    fn tc_conversions(
        &self,
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>,
    ) -> u32 {
        self.tc_conversion_pos(f, h).iter().count() as u32
    }

    /// Edits tags to included T>>C counts.
    fn push_tc_conv_aux(
        &mut self,
        auxtag: &[u8],
        f: &Option<ConvFilter>,
        h: &HashMap<u32, String>
    ) -> Result<(), NascentMolError> {
        if !self.is_possible_nascent() {
            match self.push_aux(auxtag, &Aux::Integer(0)) {
                Ok(_i) => return Ok(()),
                Err(_i) => return Err(NascentMolError::Some),
            }
        }
        let tc_aux = Aux::Integer(self.tc_conversions(f, h) as i64);
        match self.push_aux(auxtag, &tc_aux) {
            Ok(_i) => Ok(()),
            Err(_i) => Err(NascentMolError::Some),
        }
    }
}

quick_error! {
    /// Error for T>>C count annotation.
    #[derive(Debug, Clone)]
    pub enum NascentMolError {
        Some {
            description("error finding or pushing tc conv info")
        }
    }
}

/// Creates a map for converting TIDs to chromosome names.
pub fn tid_2_contig(h: &HeaderView) -> HashMap<u32, String> {
    let mut dict: HashMap<u32, String> = HashMap::with_capacity(100);
    for (i, t) in h
        .target_names()
        .iter()
        .map(|a| str::from_utf8(a).unwrap())
        .enumerate()
    {
        dict.insert(i as u32, t.to_owned());
    }
    dict
}

#[cfg(test)]
/// Tests are performed on reads mapped to MM10.
mod tests {
    use super::*;
    use handler::Nascent;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn forward_read_insertions() {
        let bampath = Path::new("test/test_insertion_forward.myc.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 2);
    }

    #[test]
    fn forward_read_deletions() {
        let bampath = Path::new("test/test_deletion_forward.myc.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 1);
    }

    #[test]
    fn forward_read_intron() {
        let bampath = Path::new("test/test_intron_forward.myc.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 5);
    }

    #[test]
    fn rev_read_insertions() {
        let bampath = Path::new("test/test_insertion_rev.brd2.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 1);
    }

    #[test]
    fn rev_read_deletions() {
        let bampath = Path::new("test/test_deletion_rev.brd2.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 2);
    }

    #[test]
    fn rev_read_intron() {
        let bampath = Path::new("test/test_intron_rev.brd2.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header().to_owned();
        let tid_lookup = tid_2_contig(&hdrv);
        let tcc: Vec<u32> = bam
            .records()
            .map(|a| a.unwrap())
            .into_iter()
            .map(|a| a.tc_conversions(&None, &tid_lookup))
            .collect();
        assert_eq!(tcc[0], 1);
    }

}
