ginger-lib: a RUST library for zk-SNARK proofs
================================================================================

Ginger-lib is a general purpose __zk-SNARK__ library that supports SNARK recursion.

Originally a fork of the [ZEXE](https://github.com/scipr-lab/ZEXE) project, ginger-lib was created with the goal of being a use-case-agnostic library, not focused on ZEXE's Decentralized Private Computation (“DPC”) or any other specific application. Ginger-lib is meant to be used specifically as a developer toolset to implement zk-SNARK schemes. The library supports full recursive proof composition and introduces several new cryptographic primitives and gadgets.

**Ginger** in Italian is "Zenzero", and this “**Zen zero**-knowledge” library was indeed developed to add some more spice to the already hot global zk-SNARK movement!


## Overview

Ginger-lib was built with the goal of developing a library that could be easily integrated and used by any project that needs to implement a zk-SNARK proving system, or a part of it. As such, it includes only zk-SNARK core objects and functionalities, with some closely related ancillary tools, and it does not assume any specific use case. 

Its first release comes with a complete set of tool to implement recursive proof composition, as detailed in the ["Scalable Zero Knowledge via Cycles of Elliptic Curves"](https://eprint.iacr.org/2014/595.pdf) paper.
In particular, it adds to the original ZEXE code the following elements:

-   __the MNT4-753 and MNT6-753 curves__
    with their Base and Scalar Fields in a full cycle. The two curves were originally sampled for and made available by the Coda project.
-   __mixed-domain FFT__
    to allow efficient conversion between coefficient and point-value polynomial representations in the domains needed to support large circuits.

The library includes also some new cryptographic primitives, implemented to be efficiently modelled in a SNARK, and in particular:

-   __"POSEIDON" hash function.__
    A SNARK-friendly hash function, i.e. well suited for being modelled as a circuit because of its low constraint need. Our [POSEIDON](https://eprint.iacr.org/2019/458.pdf) implementation has been tailored for the MNT4 and MNT6 scalar fields and heavily optimized for performance.
-   __Schnorr signature.__
    It relies on POSEIDON as random oracle, and it's fully optimized to be verified in a SNARK.
-   SNARK-friendly __Verifiable Random Function.__
    A VRF implementation based on our Schnorr and POSEIDON primitives.
-   a SNARK-friendly __Merkle Tree.__
    It uses POSEIDON as its hash function.

Ginger-lib comes also with a set of new gadgets to enforce the primitives listed above: 

-   __"proof-verifier" gadgets__ for both the MNT4 and the MNT6 curve.
    These two gadgets are the core components of the full cycle recursion. Each of them can enforce that a previous SNARK proof verifies.
-   __"POSEIDON" hash gadget.__
    A gadget enforcing that a pre-image value hashes to a defined result.
-   __Schnorr signature gadget.__
    It's a gadget that enforces that a signature, that was signed with our Schnorr primitive, verifies.
-   __VRF gadget.__
    It's a gadget that enforces that a pubblic key and a message VRF-output to a defined value.
-   __Merkle Tree gadgets.__
    A set of Merkle Tree gadgets modelling our POSEIDON-based Merkle Tree. In particular, one of the gadgets enforces that all the Merkle Tree leaves hash to a defined Merkle Root. A second gadget enforces that a value and an authentication path hash to a known Merkle Root.

Extensive automated tests have been introduced for all the added implementations.

Since it was developed to support real-world applications, ginger-lib has a strong focus on performance. Some code has already been optimized for optimal time performance; more specifically, the heaviest optimizations were performed on the implementation of the POSEIDON hash function. The number of required, expensive field inversions could be reduced by selecting "x^-1" as S-Box(x), a function compatible with the MNT4 and MNT6 curves, and then combining several of them into a "one inversion + field multiplications" operation. The same inversion trick was used also to speed up hashing in Merkle Tree processing. Further performance improvements were obtained by parallelizing the code for multi-core implementation, and by working on speeding-up the implementation of field multiplication.

Continuous performance improvement will be a key goal for all future releases and improvements of the library.  


## Directory structure

The high-level structure of the repository is as follows.

* [`algebra`](algebra): Rust crate that provides all the mathematical "bricks": finite fields, elliptic curves, FFT.
* [`primitives`](primitives): Rust crate that implements all the key cryptographic primitives.
* [`proof-systems`](proof-systems): Rust crate that implements the [Groth16](https://ia.cr/2016/260) and [GM17](https://ia.cr/2017/540) zk-SNARK proving systems.
* [`r1cs-core`](r1cs/core): Rust crate that defines core interfaces for a Rank-1 Constraint System (R1CS)..
* [`r1cs-std`](r1cs/gadgets/std): Rust crate that provides various gadgets used as building blocks of more complex R1CS circuits.
* [`r1cs-crypto`](r1cs/gadgets/crypto): Rust crate that provides various cryptographic primitives gadgets. 

In addition, there is a  [`bench-utils`](bench-utils) crate which contains infrastructure for benchmarking. It includes macros for timing code segments. It hasn't been changed from the original ZEXE implementation.

## Documentation

Detailed information about the choices made when designing and implementing our primitives and gadgets is available in the [`docs/`](docs/) directory. You can find in particular the following documents:

* [`Schnorr Verify`](docs/SchnorrVerify) documents our length-restricted variant of the Schnorr signature,  and its verification circuit
* [`VRF gadget`](docs/VRFgadget) documents our single-exponent based VRF, with Schnorr NIZK proof of correctness, and its verification circuit

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
* Please gpg sign your commits 
* Please make sure you push your pull requests to the development branch

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)

## License

ginger-lib is licensed under the following license:

 * MIT license ([LICENSE-MIT](http://opensource.org/licenses/MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in ginger-lib by you shall be licensed as above, without any additional terms or conditions.

[![License MIT](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

## Acknowledgements

The library was developed by [Horizen \(formerly ZenCash\)](https://horizen.global), as part of their effort to implement the [Zendoo](https://eprint.iacr.org/2020/123.pdf "Zendoo") sidechain system.
The project started by modifying a forked code base originally developed by the SCIPR Lab researchers for their [**ZEXE**](https://github.com/scipr-lab/ZEXE) project.  
ZEXE had previously borrowed some code from the Zcash/ECC [**Bellman**](https://github.com/zcash/librustzcash/tree/master/bellman) library.  
Some of the objects made available in this repo were adapted by the work performed by O(1) Labs for their [**Coda**](https://github.com/CodaProtocol/coda) project.  
Ginger-lib owes deeply to SCIPR Lab's [**LibSNARK**](https://github.com/scipr-lab/libSNARK), the real foundation of all practical zk-SNARK development activities. 
