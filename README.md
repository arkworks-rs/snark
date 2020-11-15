<h1 align="center">SNARK and Relation Traits</h1>

<p align="center">
    <img src="https://github.com/arkworks-rs/algebra/workflows/CI/badge.svg?branch=master">
    <a href="https://github.com/arkworks-rs/algebra/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/algebra/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/arkworks-rs/algebra"><img src="https://deps.rs/repo/github/arkworks-rs/algebra/status.svg"></a>
</p>

The arkworks ecosystem consists of Rust libraries for designing and working with __zero knowledge succinct non-interactive arguments (zkSNARKs)__. This repository contains efficient libraries that describe interfaces for zkSNARKs, as well as interfaces for programming them.

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Directory structure

This repository contains two Rust crates:

* [`ark-snark`](snark): Provides generic traits for zkSNARKs
* [`ark-relations`](relations): Provides generic traits for NP relations used in programming zkSNARKs, such as R1CS

## Overview

This repository provides the core infrastucture for using the succinct argument systems that arkworks provides. Users who want to produce arguments about various problems of interest will first reduce those problems to an NP relation, various examples of which are defined in the `ark-relations` crate. Then a SNARK system defined over that relation is used to produce a succinct argument. The `ark-snark` crate defines a `SNARK` trait that encapsulates the general functionality, as well as specific traits for various types of SNARK (those with transparent and universal setup, for instance). Different repositories within the arkworks ecosystem implement this trait for various specific SNARK constructions, such as [Groth16](https://github.com/arkworks-rs/groth16), [GM17](https://github.com/arkworks-rs/gm17), and [Marlin](https://github.com/arkworks-rs/marlin).

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the libraries:
```bash
git clone https://github.com/arkworks-rs/snark.git
cd algebra
cargo build --release
```

## Tests
This library comes with comprehensive unit and integration tests for each of the provided crates. Run the tests with:
```bash
cargo test --all
```

## License

The crates in this repo are licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in this library by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[zexe]: https://ia.cr/2018/962

## Acknowledgements

This work was supported by:
a Google Faculty Award;
the National Science Foundation;
the UC Berkeley Center for Long-Term Cybersecurity;
and donations from the Ethereum Foundation, the Interchain Foundation, and Qtum.

An earlier version of this library was developed as part of the paper *"[ZEXE: Enabling Decentralized Private Computation][zexe]"*.
