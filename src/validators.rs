use regex::{Regex};
use std::path::Path;
use std::ffi::OsStr;

pub fn tag_is_reserved_local(a: String) -> Result<(), String> {
    let local_tag_re: Regex = Regex::new(r"^[XYZ|a-z][[:alpha:]]$")
                                    .unwrap();
    match local_tag_re.is_match(&a) {
        true  => Ok(()),
        false => Err(String::from("Your tag must be a valid tag reserved for local use. See SAMTags spec sheet."))
    }
}

pub fn dir_exists(a: String) -> Result<(), String> {
    let p = Path::new(&a);

    match p.parent().unwrap().exists() {
              true => Ok(()),
              false => Err(String::from("Your intended result directory doesn't exist yet."))
           }
}

pub fn bam_seems_ok(a: String) -> Result<(), String> {
    let p = Path::new(&a);

    match p.is_file() {
        true => { match p.extension().unwrap().to_str().unwrap() == "bam" {
                true => Ok(()),
                false => Err(String::from("Input bam does not appear to be a bam."))
            }
        },
        false => Err(String::from("Input bam doesn't seem to exist."))
    }
}

pub fn vcf_seems_ok(a: String) -> Result<(), String> {
    let p = Path::new(&a);
    let ext = p.extension().unwrap().to_str().unwrap();
    let is_acceptable_ext =  ext == "vcf" || ext == "bcf";

    let tbi = format!("{}.tbi",&a);
    let cgi = format!("{}.cgi",&a);

    let idx_exists = {
        Path::new(&tbi).is_file() ||
            Path::new(&cgi).is_file()
    };

    match p.is_file() {
        true => { match is_acceptable_ext {
                true => {
                    match idx_exists {
                        true => Ok(()),
                        false => Err(String::from("Please index your vcf/bcf."))
                    }
                },
                false => Err(String::from("Input vcf/bcf does not appear to be a vcf/bcf."))
            }
        },
        false => Err(String::from("Input vcf/bcf doesn't seem to exist."))
    }
}
