use regex::{Regex};
use std::path::Path;
use std::ffi::OsStr;

/// Checks that the provided tag for annotated conversion is appropriate per SAM spec.
pub fn tag_is_reserved_local(a: String) -> Result<(), String> {
    let local_tag_re: Regex = Regex::new(r"^[XYZ|a-z][[:alpha:]]$")
                                    .unwrap();
    match local_tag_re.is_match(&a) {
        true  => Ok(()),
        false => Err(String::from("Your tag must be a valid tag reserved for local use. See SAMTags spec sheet."))
    }
}

/// Checks that the output bam path is writeable.
pub fn dir_exists(a: String) -> Result<(), String> {
    let p = Path::new(&a);

    match p.parent().unwrap().exists() {
              true => Ok(()),
              false => Err(String::from("Your intended result directory doesn't exist yet."))
           }
}

/// Checks that the input bam seems to exist.
pub fn bam_seems_ok(a: String) -> Result<(), String> {
    /// Not guaranteed to work as a validator in all cases.
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

/// Checks that input vcf/bcf seems to exist and is indexed.
pub fn vcf_seems_ok(a: String) -> Result<(), String> {
    /// Not guaranteed to work as a validator in all cases.
    let p = Path::new(&a);
    let ext = p.extension().unwrap().to_str().unwrap();
    let is_acceptable_ext =  ext == "vcf" || ext == "bcf";

    let tbi = format!("{}.tbi",&a);
    let csi = format!("{}.csi",&a);

    let idx_exists = {
        Path::new(&tbi).is_file() ||
            Path::new(&csi).is_file()
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
