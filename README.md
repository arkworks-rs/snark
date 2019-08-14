<h1 align="center">ZEXE (Zero knowledge EXEcution)</h1>

<p align="center">
    <a href="https://travis-ci.org/scipr-lab/zexe"><img src="https://travis-ci.org/scipr-lab/zexe.svg?branch=master"></a>
    <a href="https://github.com/scipr-lab/zexe/blob/master/AUTHORS"><img src="https://img.shields.io/badge/authors-SCIPR%20Lab-orange.svg"></a>
    <a href="https://github.com/scipr-lab/zexe/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
   <a href="https://github.com/scipr-lab/zexe/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>


___ZEXE___ (pronounced */zeksē/*) is a Rust library for decentralized private computation.


This library was initially developed as part of the paper *"[ZEXE: Enabling Decentralized Private Computation][zexe]"*, and it is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview

This library implements a ledger-based system that enables users to execute offline computations and subsequently produce publicly-verifiable transactions that attest to the correctness of these offline executions. The transactions contain *zero-knowledge succinct arguments* (zkSNARKs) attesting to the correctness of the offline computations, and provide strong notions of privacy and succinctness.

- **Privacy** - transactions reveal no information about the offline computation.
- **Succinctness** - transactions can be validated in time that is independent of the offline computation.
- **Application isolation** - malicious applications cannot affect the execution of honest applications.
- **Application interaction** -  applications can safely communicate with each other.

Informally, the library provides the ability to create transactions that run arbitrary (Turing-complete) scripts on hidden data stored on the ledger. In more detail, the library implements a cryptographic primitive known as *decentralized private computation* (DPC) schemes, which are described in detail in the [ZEXE paper][zexe].

## Directory structure

This repository contains several Rust crates that implement the different building blocks of ZEXE. The high-level structure of the repository is as follows.

* [`algebra`](algebra): Rust crate that provides finite fields and elliptic curves
* [`dpc`](dpc): Rust crate that implements DPC schemes (the main cryptographic primitive in this repository)
* [`r1cs-core`](r1cs-core): Rust crate that defines core interfaces for a Rank-1 Constraint System (R1CS)
* [`r1cs-std`](r1cs-std): Rust crate that provides various gadgets used to construct R1CS
* [`gm17`](gm17): Rust crate that implements the zkSNARK of [Groth and Maller][GM17]

In addition, there is a  [`bench-utils`](bench-utils) crate which contains infrastructure for benchmarking. This crate includes macros for timing code segments and is used for profiling the building blocks of ZEXE.

[GM17]: https://ia.cr/2017/540

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/scipr-lab/zexe.git
cd zexe/dpc
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test
``` 

Lastly, this library comes with benchmarks for the following crates:

- [`algebra`](algebra)
- [`dpc`](dpc)

These benchmarks require the nightly Rust toolchain; to install this, run `rustup install nightly`. Then, to run benchmarks, run the following command:
```bash
cargo +nightly bench
```

## License

ZEXE is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in ZEXE by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[zexe]: https://ia.cr/2018/962

## Reference paper

[_ZEXE: Enabling Decentralized Private Computation_][zexe]    
[Sean Bowe](https://www.github.com/ebfull), Alessandro Chiesa, Matthew Green, Ian Miers, [Pratyush Mishra](https://www.github.com/pratyush), [Howard Wu](https://www.github.com/howardwu)    
*IACR ePrint Report 2018/962*

## Acknowledgements

This work was supported by:
a Google Faculty Award;
the National Science Foundation;
the UC Berkeley Center for Long-Term Cybersecurity;
and donations from the Ethereum Foundation, the Interchain Foundation, and Qtum.

Some parts of the finite field arithmetic, elliptic curve arithmetic, FFTs, and multi-threading infrastructure in the `algebra` crate have been adapted from code in the [`ff`](https://github.com/zkcrypto/ff), [`pairing`](https://github.com/zkcrypto/pairing), and [`bellman`](https://github.com/zkcrypto/bellman) crates, developed by [Sean Bowe](https://www.github.com/ebfull) and others from Zcash.
