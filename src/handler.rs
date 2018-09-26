
use regex::{Regex, Replacer};
use rust_htslib::bam::record::Record;
use super::expander::md_expanded;

// TODO useful error if no tags, handle if unaligned


pub trait Nascent {
	fn mutated(&self) -> bool;

	fn md_ref_seq(&self) -> String;

	fn md_tag(&self) -> String;

	fn get_tc_mismatch_pos(&self) -> Vec<usize>;
}

impl Nascent for Record {
	fn mutated(&self) -> bool {
		
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

	fn get_tc_mismatch_pos(&self) -> Vec<usize> {
		

		let exp_md = self.md_ref_seq();


		vec![1,2,3]
	}
}

