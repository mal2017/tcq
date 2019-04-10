#!/bin/bash -e

rustup target add x86_64-unknown-linux-musl &&
cargo install --root $PREFIX --target x86_64-unknown-linux-musl
