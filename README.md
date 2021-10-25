<h1 align="center">ginger-lib: a RUST library for recursive SNARKs using Darlin</h1>

<p align="center">
    <a href="https://github.com/HorizenOfficial/ginger-lib/tree/master/AUTHORS"><img src="https://avatars.githubusercontent.com/u/29291571?s=20&v=4",style="width: 5vw"></a>
    <a href="https://travis-ci.com/github/HorizenOfficial/ginger-lib"><img src="https://app.travis-ci.com/HorizenOfficial/ginger-lib.svg?branch=master"></a>
    <a href="https://github.com/HorizenOfficial/ginger-lib/tree/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
   <a href="https://github.com/HorizenOfficial/ginger-lib/tree/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>


Ginger-lib is a high-performance library for building succinct zero-knowledge arguments by means of the *Darlin protocol suite*. The core piece of the protocol suite is the [Darlin](https://eprint.iacr.org/2021/930.pdf) argument system, a recursion-friendly variant of the [Marlin](https://eprint.iacr.org/2019/1047) zk-SNARK for rank-1 constraint systems (R1CS). Darlin relies on the *dlog* polynomial commitment scheme and uses an aggregation technique similar to [Halo](https://eprint.iacr.org/2019/1021) for amortizing the computational costs of both prover and verifier over time. The scheme requires no trusted setup and allows ordinary sized elliptic curves. See our reference paper for details.

The library is based on a fork from [arkworks](https://github.com/arkworks-rs/), adapted to the specific needs of the Darlin protocol suite. 

## Overview

The full protocol suite comes with a variety of proof systems to support a wide range of applications. 
These are

- *Coboundary Marlin* for simple non-recursive proofs, 
- *Darlin* for recursion, including a standard set of circuits for proof composition,
- *Rainbow Marlin*, yet another Marlin variant, which transforms Darlin proofs into ordinary Coboundary Marlin proofs.   

A detailed specification of  coboundary Marlin and Darlin, including security proofs are given in our reference paper on [Darlin](https://eprint.iacr.org/2021/930.pdf). In short, coboundary Marlin is an optimization of Marlin. It uses a simpler "sumcheck" argument, and applies a different matrix arithmetization based on the normalized Lagrange kernel. Darlin is a recursive argument (or, "accumulator SNARK") that aggregates both the dlog hard parts as well as Marlin's inner sumchecks over "time". Inner sumcheck aggregation is done cross-circuit, and a Darlin proof includes an inner sumcheck aggregator (or, *rainbow accumulator*) which supports a given pre-defined family of circuits. 

Rainbow Marlin verifies a previous Darlin proof by running a cross-circuit inner sumcheck argument for its rainbow accumulator, overall transforming Darlin proofs into ordinary Marlin proofs.

The library comes with a collection of circuits ("gadgets"), optimized for a lower R1CS density whenever needed.  

## Directory structure

The high-level structure of the repository is as follows:

* [`algebra`](algebra):  implements the arithmetics for the mathematical base components:  large integers, finite fields, elliptic curves, and fast Fourier transform.
* [`primitives`](primitives): serves basic cryptographic primitives (such as hash functions, Merkle trees, signature schemes, verifiable random functions).
* [`proof-systems`](proof-systems): The main crate for the Darlin protocol suite. Provides the traits and structs for proof carrying data and the above mentioned proof systems. 
* [`r1cs-core`](r1cs/core): Defines core interfaces for rank-1 constraint systems.
* [`r1cs-std`](r1cs/gadgets/std): This crate provides  elementary "standard" circuits ("gadgets"): Boolean operations, native field and elliptic curve arithmetics.  
* [`r1cs-crypto`](r1cs/gadgets/crypto): provides the circuits for various cryptographic primitives, such as the Poseidon hash, signature schemes, and SNARK verifiers.

In addition, there is a  [`bench-utils`](bench-utils) crate which contains infrastructure for benchmarking, including macros for timing code segments. 

## Release Note

The current release serves does not yet provide proof composition. The [`proof-systems`](proof-systems) subcrate [`darlin`](proof-systems/src/darlin) prepares for the full Darlin protocol suite by providing the traits and structs necessary for proof carrying data, and puts simple Coboundary Marlin proofs into this framework. It further provides additional tools around the verification of future Darlin along ordinary Marlin proofs at large scale:
- A batch verifier for Darlin/Marlin proofs, and 
- a post-processor for batches of Darlin/Marlin proofs, which aggregates their dlog hard parts into a single one.

## Build instructions

The library compiles on the `1.51.0 stable` toolchain of the Rust compiler. 
```bash
git clone https://github.com/HorizenOfficial/ginger-lib.git
cd ginger-lib
cargo build --release
```
Run tests using
```bash
cargo test --all-features 
```
More detailed information can be found in our [build guide](build_guide.md).


## License

This library is licensed under either of the following licenses, at your discretion.

 * [Apache License Version 2.0](LICENSE-APACHE)
 * [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.
