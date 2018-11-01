use regex::{Regex};

pub fn tag_is_reserved_local(a: String) -> Result<(), String> {
    let local_tag_re: Regex = Regex::new(r"^[XYZ|a-z][[:alpha:]]$")
                                    .unwrap();

    match local_tag_re.is_match(&a) {
        true  => Ok(()),
        false => Err(String::from("Your tag must be a valid tag reserved for local use. See SAMTags spec sheet."))
    }
}
