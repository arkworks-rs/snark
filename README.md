<h1 align="center">SNARK and Relation Traits</h1>

<p align="center">
    <img src="https://github.com/arkworks-rs/algebra/workflows/CI/badge.svg?branch=master">
    <a href="https://github.com/arkworks-rs/algebra/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/algebra/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/arkworks-rs/algebra"><img src="https://deps.rs/repo/github/arkworks-rs/algebra/status.svg"></a>
</p>

The **arkworks** ecosystem consists of Rust libraries for designing and working with **zero-knowledge succinct non-interactive arguments (zkSNARKs)**. This repository provides efficient, modular traits for zkSNARKs and **NP relations**, forming the foundation for building and verifying zero-knowledge proofs.

This library is released under the **MIT License** and the **Apache v2 License** (see [License](#license)).

> 🚨 **WARNING**: This is an **academic proof-of-concept** and has **not been audited** for security. Do **not** use it in production environments that handle sensitive data.

---


## Directory Structure
This repository contains two Rust crates:

- **`ark-snark`** – Defines **generic traits** for zkSNARKs.
- **`ark-relations`** – Provides **NP relation traits** used in programming zkSNARKs, such as **R1CS**.

---
## Overview

This repository provides the core infrastructure for using the succinct argument systems that arkworks provides. Users who want to generate proofs for different computational problems must first reduce them to an NP relation. Examples of such relations are provided in the `ark-relations` crate. Then a SNARK system defined over that relation is used to produce a succinct argument. The `ark-snark` crate defines a `zkSNARK` trait that encapsulates the general functionality, as well as specific traits for various types of SNARK (those with transparent and universal setup, for instance). Other repositories in the arkworks ecosystem implement this trait for specific SNARK constructions, such as [Groth16](https://github.com/arkworks-rs/groth16), [GM17](https://github.com/arkworks-rs/gm17), and [Marlin](https://github.com/arkworks-rs/marlin).

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

Use `cargo`, the standard Rust build tool, to build the libraries:
```bash
git clone https://github.com/arkworks-rs/snark.git
cd snark
cargo build --release
```

## Testing
This library comes with comprehensive unit and integration tests for each of the provided crates. To run the tests, use:
```bash
cargo test --all
```
## Related Libraries
This repository is part of the **arkworks ecosystem**. Other relevant libraries include:

- [`ark-crypto-primitives`](https://github.com/arkworks-rs/crypto-primitives) – Cryptographic primitives for zkSNARKs.
- [`ark-r1cs-std`](https://github.com/arkworks-rs/r1cs-std) – Standardized constraints for R1CS.
- [`ark-groth16`](https://github.com/arkworks-rs/groth16) – Groth16 zkSNARK implementation.
- [`ark-marlin`](https://github.com/arkworks-rs/marlin) – Marlin zkSNARK implementation.

These libraries **extend the functionality** of `ark-snark` and `ark-relations`, making it easier to construct and verify zkSNARK proofs.

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
