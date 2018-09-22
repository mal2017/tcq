use regex::Regex;
use rust_htslib::bam::record::Record;

// TODO useful error if no tags, handle if unaligned

pub fn has_var(md: &String) -> bool {
	let nums_re =  Regex::new(r"\D").unwrap();
	nums_re.is_match(md)
	// circle seq, insert biotinylated bases
	// circularize and ligate ends
	// chew up with enzyme from circle seq
}

pub trait Nascent<T> {
	fn mutated(&self) -> bool;

	fn ref_seq(&self) -> String;

	fn md_tag(&self) -> String;
}

impl Nascent<Record> for Record {
	fn mutated(&self) -> bool {
		
		self.aux(b"NM").unwrap().integer() > 0
	}

	fn md_tag(&self) -> String {
		String::from_utf8(self.aux(b"MD")
							  .unwrap()
							  .string()
							  .to_owned())
							  .unwrap()
	}

	fn ref_seq(&self) -> String {
		self.md_tag()
		// get idx
	}

}


//fn md2ref() {
//	let re = Regex::new()
//}