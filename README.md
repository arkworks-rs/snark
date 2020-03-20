ginger-lib: a RUST library for zk-SNARK proofs
================================================================================

This library implements __zk-SNARK__ schemes, which are a cryptographic method
for proving/verifying, in zero knowledge, the integrity of computations.

The library was initially developed and is maintained by the Zen Blockchain Foundation, the legal entity of [Horizen \(formerly ZenCash\)](https://horizen.global), as part of their implementation effort of the [Zendoo](https://eprint.iacr.org/2020/123.pdf "Zendoo") sidechain system. The code base started as a fork of the [Zexe](https://github.com/scipr-lab/zexe) project. It is released under the Apache and MIT License, and offered to the developer community as a contribution to the larger adoption and the improvement of SNARG technologies, in and beyond the crypto world.

Ginger in Italian is "Zenzero". And **Zen** has indeed developed this **zero**-knowledge library to further spice-up global zk-SNARK development!  


## Overview

Ginger-lib was built with the goal of developing a library that could be easily integrated and used by any project that needs to implement a zk-SNARK proving system, or a part of it. As such, it includes only zk-SNARK core objects and functionalities, with some closely related ancillary tools, and it does not assume any specific use case. 

Its first release comes with a complete set of tool to implement recursive proof composition, as detailed in the ["Scalable Zero Knowledge via Cycles of Elliptic Curves"](https://eprint.iacr.org/2014/595.pdf) paper.
In particular, it adds to the original Zexe code the following elements:

-   __the MNT4-753 and MNT6-753 curves__
    with their Base and Scalar Fields in a full cycle
-   __mixed-domain FFT__
    to allow efficient conversion between coefficient and point-value polynomial representations in the domains needed to support large circuits

The library includes also some new cryptographic primitives, implemented to be efficiently modelled in a Snark, and in particular:

-   __"Poseidon" hash function.__
    A snark-friendly, ideal to be modelled in a circuit for its low constraint need. Our [Poseidon](https://eprint.iacr.org/2019/458.pdf) implementation has been heavily optimized for performance. 
-   __Schnorr signature.__
    It relies on Poseidon as random oracle, and it's fully optimized to be verified in a Snark.
-   snark-friendly __Verifiable Random Function.__
    A VRF based on our Schnorr and Poseidon implementations.
-   a snark-friendly __Merkle Tree.__
    It uses Poseidon as its hash function.

Ginger-lib comes also with a list of new gadgets that can model the above listed primitives: 

-   __"proof-verifies" gadget for both the MNT4 and the MNT6 curve.__
    These two circuits are the core components of full cycle recursion. They can snark-prove that a previous snark proof verifies.
-   __"Poseidon" hash gadget.__
    A circuit that can prove that a pre-image witness hashes to a known public input value.
-   __Schnorr signature gadget.__
    It's a circuit that can snark-prove that a known signature, that was signed with our Schnorr primitive, verifies.
-   __VRF gadget.__
    It's a circuit that can prove that a pubkey and a messages, as witnesses, output to a known VRF output. This VRF output is a public input.
-   __Merkle Tree gadgets.__
    Merkle Tree circuits adapted and optimised for our Poseidon-based Merkle Tree. In particular, there is a gadget that prove that all the leaves (witnesses) hash to a known Merkle Root (public input), and a circuit that can prove that a value and an authentication path (witnesses) hash to a known Merkle Root (public input).

Also, extensive automated tests have been introduced for all the added implementations.

Since it was developed to support real-world applications, ginger-lib has a strong focus on performance; some code has already been optimized for optimal time performance, and continuous performance improvement will be a key goal of the core development team. 

The heaviest optimizations were performed on the Poseidon implementation; by selecting "x^-1" as S-Box(x), a function compatible with the MNT4 and MNT6 curves, the number of required expensive field inversions could be reduced by processing them in bulks. The same inversion trick was used also to speed up hashing in Merkle Tree processing. Further performance improvements were obtained by parallelization for multi-core implementation, and by working on the implementation of the field multiplication.
 


## Directory structure

The high-level structure of the repository is as follows.

* [`algebra`](algebra): Rust crate that provides all the mathematical "bricks": finite fields, elliptic curves, fft.
* [`primitives`](primitives): Rust crate that implements all the key cryptographic primitives.
* [`proof-systems`](proof-systems): Rust crate that implements the [Groth16](https://ia.cr/2016/260) and [GM17](https://ia.cr/2017/540) zk-SNARK proving systems.
* [`r1cs-core`](r1cs/core): Rust crate that defines core interfaces for a Rank-1 Constraint System (R1CS)
* [`r1cs-std`](r1cs/gadgets/std): Rust crate that provides various gadgets used as building blocks of more complex R1CS circuits
* [`r1cs-crypto`](r1cs/gadgets/crypto): Rust crate that provides various cryptographic primitives gadgets 

In addition, there is a  [`bench-utils`](bench-utils) crate which contains infrastructure for benchmarking. This crate includes macros for timing code segments, and it hasn't been changed from the original Zexe implementation.


## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/.../ginger-lib.git
cd ginger-lib
cargo build 
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test --all-features 
``` 

Lastly, this library comes with benchmarks for the following crates:

- [`algebra`](algebra)

These benchmarks require the nightly Rust toolchain; to install this, run `rustup install nightly`. Then, to run benchmarks, run the following command:
```bash
cargo +nightly bench --all-features 
```

## Contributing



Contributions are welcomed! Bug fixes and new features can be initiated through GitHub pull requests. To speed the code review process, please adhere to the following guidelines:

* Follow code_of_conduct.md
* Follow the styling guide at doc/developer-notes.md
* Please gpg sign your commits, not required but a nice to have
* Please make sure you push your pull requests to the development branch

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)

## License

ginger-lib is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in ginger-lib by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[![License APACHE]("https://img.shields.io/badge/license-APACHE-blue.svg")](http://www.apache.org/licenses/LICENSE-2.0)
[![License MIT]("https://img.shields.io/badge/license-MIT-blue.svg")](http://opensource.org/licenses/MIT)

## Acknowledgements

This work is supported by the **Zen Blockchain Foundation**. 
The project started by modifying a forked code base originally developed by the SCIPR Lab researcher for their [**Zexe**](https://github.com/scipr-lab/zexe) project. Zexe had previously borrowed some code from Zcash/ECC [**Bellman**](https://github.com/zcash/librustzcash/tree/master/bellman) library.
Some of the objects made available in this repo were adapted by the work performed by O(1) Labs for their [**Coda**](https://github.com/CodaProtocol/coda) project.
Ginger-lib, as any of the above project, owes deeply to SCIPR Lab's [**Libsnark**](https://github.com/scipr-lab/libsnark), the real foundation of all practical zk-SNARK development. 
