
use regex::{Regex, Replacer};
use rust_htslib::bam::record::{Record, Seq};
use super::expander::md_expanded;

// TODO useful error if no tags, handle if unaligned

pub trait Nascent {
	fn is_possible_nascent(&self) -> bool;

	fn md_ref_seq(&self) -> String;

	fn md_tag(&self) -> String;

	fn cand_tc_mismatch_ref_pos(&self) -> Vec<u32>;

	fn cand_tc_mismatch_pos_tuples(&self) -> Vec<(u32,u32)>;

	fn tc_conversion_pos(&self) -> Vec<((u32,u32),bool)>;

	fn tc_conversions(&self) -> u32;
}

impl Nascent for Record {
	fn is_possible_nascent(&self) -> bool {
		
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

	fn tc_conversion_pos(&self) -> Vec<((u32,u32),bool)> {
		let cand_pos_tuples = self.cand_tc_mismatch_pos_tuples();

		let read_seq = self.seq();

		let conv_target: u8 = if self.is_reverse() {
			b'G'
		} else {
			b'C'
		};

		let enc_base_hit_itr = cand_pos_tuples.iter()
												.map(|a| a.0 as usize) // get read pos
												.map(|a| read_seq.encoded_base(a))
												.map(|a| a == conv_target);

		cand_pos_tuples.clone().into_iter()
						.zip(enc_base_hit_itr)
						.filter(|a| a.1)
						.collect()

	}

	fn tc_conversions(&self) -> u32 {
		self.tc_conversion_pos().iter()
								.filter(|a| a.1)
								.count() as u32
	}
}

