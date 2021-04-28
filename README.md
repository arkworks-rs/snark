ginger-lib: a RUST library for zk-SNARKs
================================================================================

Ginger-lib is a general purpose __zk-SNARK__ library that supports recursive composition of [Groth16](https://eprint.iacr.org/2016/260.pdf) arguments.

Originally a fork of the [ZEXE](https://github.com/scipr-lab/ZEXE) project, ginger-lib was created with the goal of being an independent library, i.e. not linked to the Decentralized Private Computation (“DPC”) application or any other specific use case. Designed as a developer toolset for implementing zero-knowledge SNARKs, ginger-lib comes with an increased collection of cryptographic primitives and matching *gadgets*, that can serve as building blocks for application-tailored statements/circuits. 

**Ginger** in Italian is "Zenzero", and this “**Zen zero**-knowledge” library was indeed developed to add some more spice to the already hot global zk-SNARK movement!


## Overview

Ginger-lib was built with the goal of being easily integrated and used by any project that needs to implement its own, application-tailored zk-SNARK. As such, it provides SNARK core objects and functionalities, and a few closely related ancillary tools. 

In particular, its first release comes with a complete set of tools for recursive proof composition as in Ben-Sasson, et al. ["Scalable Zero Knowledge via Cycles of Elliptic Curves (2014)"](https://eprint.iacr.org/2014/595.pdf). Specifically, it adds to the original ZEXE code the following elements:

-   __MNT4-753 and MNT6-753 curves__, 
    a re-implementation of [Coda](https://coinlist.co/build/coda/)'s MNT4-MNT6 cycle of pairing friendly elliptic curves for a security level of 128 bit. All curve parameters were re-checked, the pairing engine ported to Rust, the gadget collection extended with all needed recursive argument evaluation components.
-   __mixed-domain FFT__,
    to allow efficient conversion between coefficient and point-value polynomial representations in the domains needed to support large circuits.

The library includes also some additional cryptographic primitives, implemented to be efficiently modelled in a SNARK, and in particular:

-   the __POSEIDON hash function__ - 
    thanks to its efficient description as an arithmetic circuit, the [POSEIDON](https://eprint.iacr.org/2019/458.pdf) hash family is ideal for SNARKs. Our implementations for both the MNT4-753 and MNT6-753 scalar fields use the modular inversion S-Box, apply a security level of 128 bits and are optimized for performance.
-   __Schnorr NIZK proof and signature scheme__ - 
    Schnorr-like non-interactive zero-knowledge (NIZK) proof and the Schnorr signature scheme, using POSEIDON as random oracle and adapted to be efficiently integrated in a SNARK.
-   a SNARK-friendly __Verifiable Random Function (VRF)__, 
    based on our Schnorr and POSEIDON primitives.
-   a SNARK-friendly __Merkle Tree__,
    using POSEIDON as its hash function.

Along with the above primitives, ginger-lib comes with the following new gadgets: 

-   __Groth16 verification gadgets__ for both the MNT4 and the MNT6 curve.
    These  gadgets are the core components of recursive proof evaluation. They enforce that a Groth16 SNARK, based on one of these two curves, verifies.
-   __POSEIDON hash gadget__, 
    enforcing that some pre-image hashes to a given fingerprint.
-   __Schnorr proof / signature verification gadgets__,
    enforcing that a single-exponent Schnorr NIZK proof or Schnorr signature, that was created by our corresponding primitives, verifies.
-   __VRF verification gadget__,
    enforcing that a public key and message, and a VRF-output, are consistent.
-   __Merkle Tree gadget__, 
    enforcing that the authentication path of a leaf is consistent with a given Merkle root.

Extensive automated tests have been introduced for the added implementations.

Since it was developed to support real-world applications, ginger-lib has a strong focus on performance. Some of the code is already optimized in timing. More specifically, optimizations were performed on the implementation of the POSEIDON hash function:  batch hashing/verification was significantly sped up by using a "single inversion + field multiplications" replacement for multiple parallel inversions. The same trick was also used to speed up hashing in Merkle trees. Further performance improvements were obtained by parallelizing the code for multi-core implementation, and by speeding up the implementation of field multiplication.

Continuous performance improvement will be a key goal for all future releases and improvements of the library.  


 
**Please note: the library is in development, its core code is still being modified, and it has not yet undergone extensive review or testing. No guarantees are provided about its security and functionality. At this stage, it is not suitable for use in production systems.**

## Directory structure

The high-level structure of the repository is as follows:

* [`algebra`](algebra): Rust crate that provides all the mathematical "bricks": finite fields, elliptic curves, FFT.
* [`primitives`](primitives): Rust crate that implements all the key cryptographic primitives.
* [`proof-systems`](proof-systems): Rust crate that implements the [Groth16](https://ia.cr/2016/260) and [GM17](https://ia.cr/2017/540) zk-SNARK proving systems.
* [`r1cs-core`](r1cs/core): Rust crate that defines core interfaces for a Rank-1 Constraint System (R1CS).
* [`r1cs-std`](r1cs/gadgets/std): Rust crate that provides various gadgets used as building blocks of more complex R1CS circuits.
* [`r1cs-crypto`](r1cs/gadgets/crypto): Rust crate that provides various cryptographic primitives gadgets. 

In addition, there is a  [`bench-utils`](bench-utils) crate which contains infrastructure for benchmarking. It includes macros for timing code segments. It hasn't been changed from the original ZEXE implementation.

## Documentation

Detailed information about the choices made when designing and implementing our primitives and gadgets is available in the [`doc/`](doc/) directory. In particular you can find the following documents:

* [`PoseidonAndGadgets`](doc/Poseidon.pdf), it documents the parameters for our POSEIDON hash function and its verification circuit.
* [`SchnorrAndGadgets`](doc/SchnorrSignature.pdf), it explains our length-restricted variant of the Schnorr signature, and its verification circuit.
* [`SchnorrVerdictGadget`](doc/SchnorrVerdict.pdf), it describes a slight generalization of the Schnorr verification gadget, a circuit which enforces a boolean input (the "verdict") to encode the validity/non-validity of a given Schnorr signature.


## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/.../ginger-lib.git
cd ginger-lib
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test --all-features 
``` 

By default, ```cargo test``` will execute the tests concurrently on all available cores. Since some tests are resource-intensive, this may abort the tests execution. If this happens, you may want to reduce the number of cores running the tests with the command:

```bash
cargo test --all-features -- --test-threads=#threads
``` 
 Also, some tests take a long time to be executed. The reason for this is the dev profile set in the root Cargo.toml:
 
 ```
[profile.dev]
opt-level = 0
``` 

Set it to 3 for maximum speed but longer compilation times: this is suggested for executing all the tests in the project, but for single test's execution might be unnecessary.

This library comes with benchmarks for the [`algebra`](algebra) crate.
These benchmarks require the nightly Rust toolchain; to install this, run `rustup install nightly`. Then, to run benchmarks, run the following command: 
```bash
cargo +nightly bench --all-features 
```

Compiling with `adcxq`, `adoxq` and `mulxq` instructions can lead to a 30-70% speedup. These are available on most `x86_64` platforms (Broadwell onwards for Intel and Ryzen onwards for AMD). Run the following command:
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly test/build/bench --features llvm_asm
```
Tip: If optimising for performance, your mileage may vary with passing `--emit=asm` to `RUSTFLAGS`.

To bench `algebra-benches` with greater accuracy, especially for functions with execution times on the order of nanoseconds, use the `n_fold` feature to run selected functions 1000x per iteration. To run with multiple features, make sure to double quote the features.
```bash
cargo +nightly bench --features "n_fold"
```
__Note:__ Some of the dependencies between the crates in GingerLib are specified via Git rather than via local paths: this is due to a cross-dependency issue between GingerLib's crates and some external crates. 
One example of such errors is in crate `proof-systems`: it depends both on `algebra` and on external crates located in [marlin](https://github.com/HorizenLabs/marlin) and [poly-commit](https://github.com/HorizenLabs/poly-commit) depending on `algebra` too; if the version of `algebra` on which these crates depend is not exactly the same, a compilation error will occur:
```bash
error[E0308]: mismatched types [...] note: perhaps two different versions of crate `algebra` are being used?
```
By specifying in all the crates the dependency on `algebra` in Git form, we ensure that all the crates will take the same version; however, if during development `algebra` crate is modified, we would be forced to push the changes to Git first before seeing them applied in local. For this reason, in the root `Cargo.toml`, we pushed instructions allowing to override Git dependencies with (local) path dependencies; unfortunately, this will require to store locally all the crates involved in the cross-dependency issue and to  comment/uncomment these lines (if needed) before/after pushing changes.
We are considering to restructure the involved repositories to avoid this issue.

## Contributing

Contributions are welcomed! Bug fixes and new features can be initiated through GitHub pull requests. To speed the code review process, please adhere to the following guidelines:

* Follow Horizen repositories' *code of conduct*
* Follow Horizen repositories' *styling guide* 
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
The project started by modifying a forked code-base originally developed by the SCIPR Lab researchers for their [**ZEXE**](https://github.com/scipr-lab/ZEXE) project.  
ZEXE had previously borrowed some code from the Zcash/ECC [**Bellman**](https://github.com/zcash/librustzcash/tree/master/bellman) library.  
Some of the objects made available in this repo were adapted by the work performed by O(1) Labs for their [**Coda**](https://github.com/CodaProtocol/coda) project.  
Ginger-lib owes deeply to SCIPR Lab's [**LibSNARK**](https://github.com/scipr-lab/libSNARK), the real foundation of all practical zk-SNARK development activities. 
