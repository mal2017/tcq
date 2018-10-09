
use regex::{Regex};
use rust_htslib::bam::record::{Record, Aux};
use super::expander::md_expanded;
use super::filter::ConvFilter;
use std::collections::HashMap;
use std::str;
use rust_htslib::bam::HeaderView;
use std::ops::Range;
use bio::utils::Interval;

pub trait Nascent {
	fn is_possible_nascent(&self) -> bool;

	fn md_ref_seq(&self) -> String;

	fn md_tag(&self) -> String;

	fn cand_tc_mismatch_ref_pos(&self) -> Vec<u32>;

	fn cand_tc_mismatch_pos_tuples(&self) -> Vec<(u32,u32)>;

	fn tc_conversion_pos(&self, f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> Vec<((u32,u32),bool)>;

	fn tc_conversions(&self, f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> u32;

	fn push_tc_conv_aux(&mut self, auxtag: &[u8], f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> Result<(), NascentMolError>;
}

impl Nascent for Record {
	fn is_possible_nascent(&self) -> bool {
		// TODO check reads have md tags,unless unmapped
		// TODO check bits for revcomp set on - strand reads
		if self.is_unmapped() { return false };

		lazy_static! { // for speeeeed
			static ref revr_re: Regex  = Regex::new(r"A").unwrap();
			static ref forw_re: Regex  = Regex::new(r"T").unwrap();
		}

		let md = self.md_tag();

		match self.is_reverse() {
			true => {
				revr_re.is_match(&md) //A>>G
			},
			false => {
				forw_re.is_match(&md) //T>>C
			},
		}
	}

	fn md_tag(&self) -> String {
		String::from_utf8(self.aux(b"MD")
							  .unwrap()
							  .string()
							  .to_owned())
							  .unwrap()
	}

	fn md_ref_seq(&self) -> String {
		md_expanded(self.md_tag())
	}

	fn cand_tc_mismatch_ref_pos(&self) -> Vec<u32> {
		lazy_static! { // for speeeeed
			static ref a_re: Regex  = Regex::new(r"A").unwrap();
			static ref t_re: Regex  = Regex::new(r"T").unwrap();
		}

		let exp_md = self.md_ref_seq();

		match self.is_reverse() {
			true => {
				a_re.find_iter(&exp_md).map(|a| a.start() as u32).collect() //A>>G
			},
			false => {
				t_re.find_iter(&exp_md).map(|a| a.start() as u32).collect() //T>>C
			},
		}

	}

	fn cand_tc_mismatch_pos_tuples(&self) -> Vec<(u32,u32)> {
		let cig = self.cigar();
		let start = self.pos() as u32;
		let cand_ref_pos = self.cand_tc_mismatch_ref_pos();

		cand_ref_pos.iter()
					.map(|a| cig.read_pos(start + a,true,true))
					.map(|a| a.unwrap().unwrap())
					.zip(cand_ref_pos.clone().into_iter())
					.collect()
	}

	fn tc_conversion_pos(&self, f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> Vec<((u32,u32),bool)> {
		// TODO vcf/bcf filtering here
		// TODO check oncoordinates for dels, ins, am i doing it right when i retu
		// retrieve read seq from cigar
		let mut cand_pos_tuples = self.cand_tc_mismatch_pos_tuples();
		let mut rng: Range<u32> = Range {start:0,end:0};
		let mut chr = h.get(&(self.tid() as u32)).unwrap();
		let it;
		cand_pos_tuples = match f {
			Some(k) => {
				it = &k.inner;
				//println!("{:?}", it.keys());
				cand_pos_tuples.into_iter()
							   .filter(|a| {
								   rng = Range { start: a.1, end: a.1 + 1 };
								   match it.get(chr) {
									   Some(j) => {j.find(&rng).any(|a| true)},
									   None => true,
								   }
								   // filter for blacklist overlap
							   })
								.collect()
			},
			None => cand_pos_tuples,
		};

		let read_seq = self.seq();

		let conv_target: u8 = if self.is_reverse() {
			b'G'
		} else {
			b'C'
		};

		let enc_base_hit_itr = cand_pos_tuples.iter()
												.map(|a| a.0 as usize) // get read pos
												.map(|a| read_seq.as_bytes()[a])
												.map(|a| a == conv_target);

		cand_pos_tuples.clone().into_iter()
						.zip(enc_base_hit_itr)
						.filter(|a| a.1)
						.collect()

	}

	fn tc_conversions(&self, f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> u32 {
		self.tc_conversion_pos(f, h).iter()
								.count() as u32
	}


	fn push_tc_conv_aux(&mut self, auxtag: &[u8],f: &Option<ConvFilter>, h: &HashMap<u32, String>) -> Result<(), NascentMolError> {

		if !self.is_possible_nascent() {
			match self.push_aux(auxtag, &Aux::Integer(0)) {
				Ok(_i) => {
					return Ok(())
				},
				Err(_i) => {
					return Err(NascentMolError::Some)
				},
			}
		}
		let tc_aux = Aux::Integer(self.tc_conversions(f, h) as i64);
		match self.push_aux(auxtag, &tc_aux) {
			Ok(_i) => {
				Ok(())
			},
			Err(_i) => {
				Err(NascentMolError::Some)
			},
		}
	}

}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum NascentMolError {
        Some {
            description("error finding or pushing tc conv info")
        }
    }
}


pub fn tid_2_contig(h: &HeaderView) -> HashMap<u32, String> {
	let mut dict: HashMap<u32, String> = HashMap::with_capacity(100);
	for (i,t) in h.target_names()
				  .iter().map(|a| str::from_utf8(a).unwrap())
				  .enumerate() {
		dict.insert(i as u32, t.to_owned());
	}
	dict
}
