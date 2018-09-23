use regex::{Regex, Replacer};
use regex::RegexSet;
use rust_htslib::bam::record::Record;

// TODO useful error if no tags, handle if unaligned

pub fn has_var(md: &String) -> bool {
	let nums_re =  Regex::new(r"\D").unwrap();
	nums_re.is_match(md)
	// circle seq, insert biotinylated bases
	// circularize and ligate ends
	// chew up with enzyme from circle seq
}

pub trait Nascent {
	fn mutated(&self) -> bool;

	fn ref_seq(&self) -> String;

	fn md_tag(&self) -> String;
}

impl Nascent for Record {
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
		md_expanded(&mut self.md_tag());
		// get idx
		"a".to_string()
	}

}

fn md_expanded(md: &String) {

	replace_matches(md);
}


fn replace_matches(md: &String) {
	let mut orig = md.clone();
	lazy_static! { // for speeeeed
		static ref re_perf: Regex  = Regex::new(r"\d+").unwrap();
	}
	let mut m;
	let mut expanded = false;
	let mut n: usize;
	let mut exp;
	let mut new;

	while !expanded  {
		m = re_perf.find(&orig).unwrap();
		n = m.as_str().parse().unwrap();
		exp = "M".repeat(n);
		new = orig.clone();
		replace_match1(&mut orig,m.start(),m.end(),&exp);
		//orig = new;
		//println!("{:?} {:?}",md,orig);
		expanded = true;
	}
	
}

fn replace_match1(orig: &mut String, start: usize, end: usize, repl: &String) {
	orig.replace_range(start..end,repl);
}