extern crate bio;
extern crate regex;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate log;
extern crate core;
extern crate env_logger;

/// Runs the executable functions of `tcq`.
pub mod runner;

/// Hold the interface to the engine behind `tcq`.
pub mod handler;

/// Utils for reconstructing pseudo-references from MD tags.
pub mod expander;

/// Utils for filtering RNA-seq reads.
pub mod filter;

/// Functions for validating arguments to the `tcq` executable.
pub mod validators;

/// Hacked functions for handling spliced reads.
pub mod spliced_read_utils;
