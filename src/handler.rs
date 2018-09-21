use regex::Regex;


pub fn has_var(md: &String) -> bool {
	let nums_re =  Regex::new(r"\D").unwrap();
	nums_re.is_match(md)
	// circle seq, insert biotinylated bases
	// circularize and ligate ends
	// chew up with enzyme from circle seq


}


pub fn nt2c() {

}