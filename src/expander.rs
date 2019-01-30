use regex::Regex;

/// From an MD tag reconstruct a pseudo-reference.
pub fn md_expanded(md: String) -> String {
    info!("{:?}", replace_matches(replace_dels(md.clone())));
    replace_matches(replace_dels(md))
}

/// Replace perfect matches with `M` character within a CIGAR string.
fn replace_matches(md: String) -> String {
    let orig = md.clone();
    lazy_static! { // for speeeeed
        static ref re_perf: Regex  = Regex::new(r"\d+").unwrap();
    }
    let m;
    let n: usize;
    let exp;
    let new: String; // = orig.clone();

    m = re_perf.find(&md);

    match m {
        Some(i) => {
            n = i.as_str().parse().unwrap();
            exp = "M".repeat(n);
            new = replaced_match1(&md, i.start(), i.end(), &exp);
            replace_matches(new)
        }
        None => orig,
    }
}

/// Expand an individual stretch of perfect matches.
fn replaced_match1(orig: &String, start: usize, end: usize, repl: &String) -> String {
    // TODO try not cloning here
    let mut new = orig.clone();
    new.replace_range(start..end, repl);
    new
}

/// Replace deletions with `D` characters.
fn replace_dels(md: String) -> String {
    let orig = md.clone();
    lazy_static! { // for speeeeed
        static ref re_del: Regex  = Regex::new(r"\^\D+").unwrap();
    }
    let m;
    let n: usize;
    let exp;
    let new: String; // = orig.clone();

    m = re_del.find(&md);

    match m {
        Some(i) => {
            n = i.as_str().to_owned().len();
            exp = "D".repeat(n - 1);
            new = replaced_del1(&md, i.start(), i.end(), &exp);
            replace_dels(new)
        }
        None => orig,
    }
}

/// Expand a stretch of deletions.
fn replaced_del1(orig: &String, start: usize, end: usize, repl: &String) -> String {
    // TODO try not cloning here
    let mut new = orig.clone();
    new.replace_range(start..end, repl);
    new
}
