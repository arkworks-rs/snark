# ginger-lib: a RUST library for recursive SNARKs using Darlin

<p align="center">
    <a href= "https://github.com/HorizenOfficial/ginger-lib/releases"><img src="https://img.shields.io/github/release/HorizenOfficial/ginger-lib.svg"></a>
    <a href="AUTHORS"><img src="https://img.shields.io/github/contributors/HorizenOfficial/ginger-lib.svg?"></a>
    <a href="https://travis-ci.com/github/HorizenOfficial/ginger-lib"><img src="https://app.travis-ci.com/HorizenOfficial/ginger-lib.svg?branch=master"></a>
    <a href="LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="CONTRIBUTING.md"><img src="https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square"></a>
</p>


Ginger-lib (in Italian *zen-zero*) is a high-performance library for building succinct zero-knowledge arguments by means of the *Darlin protocol suite*. The core piece of the protocol suite is the [Darlin](https://eprint.iacr.org/2021/930) argument system, a recursion-friendly variant of the [Marlin](https://eprint.iacr.org/2019/1047) zk-SNARK for rank-1 constraint systems (R1CS). Darlin relies on the *dlog* polynomial commitment scheme and uses an aggregation technique similar to [Halo](https://eprint.iacr.org/2019/1021) for amortizing the computational costs of both prover and verifier over time. The scheme requires no trusted setup and allows ordinary sized elliptic curves. See our reference paper [HGB](https://eprint.iacr.org/2021/930) for details.

The library is based on a fork from [arkworks](https://github.com/arkworks-rs/), and is adapted to the specific needs of the Darlin protocol suite. 

## Overview

The full protocol suite comes with a variety of proof systems to support a wide range of applications. 
These are

- *Coboundary Marlin* for simple non-recursive proofs, 
- *Darlin* for recursion, including a standard set of circuits for proof composition,
- *Rainbow Marlin*, yet another Marlin variant, which transforms Darlin proofs into ordinary Coboundary Marlin proofs.   

A detailed specification of  Coboundary Marlin and Darlin, including security proofs are given in [HGB](https://eprint.iacr.org/2021/930.pdf). In short, Coboundary Marlin is an optimization of Marlin. It uses a simpler "sumcheck" argument, and applies a different matrix arithmetization based on the normalized Lagrange kernel. Darlin is Coboundary Marlin turned into a recursive argument (or, *accumulator SNARK*) which aggregates both the dlog hard parts as well as Marlin's inner sumchecks over "time". Inner sumcheck aggregation is done across circuits, and a Darlin proof includes an inner sumcheck aggregator (or, *Rainbow Accumulator*) which supports a given pre-defined family of circuits. 

Rainbow Marlin is used to verify a previous Darlin proof by running a cross-circuit inner sumcheck argument for its Rainbow Accumulator, overall transforming Darlin proofs into simple Marlin proofs.

The library comes with a collection of circuits, manually optimized for a lower R1CS density whenever needed.  

## Directory structure

The high-level structure of the repository is as follows:

* [`algebra`](algebra):  implements the mathematical base components:  large integers, finite fields, elliptic curves, and fast Fourier transform.
* [`primitives`](primitives): serves basic cryptographic primitives (such as hash functions and Merkle trees, signature schemes, verifiable random functions).
* [`proof-systems`](proof-systems): This is the main crate for the Darlin protocol suite. It provides the traits and structs for proof carrying data and the above mentioned proof systems. [Groth16](https://ia.cr/2016/260) and [GM17](https://ia.cr/2017/540) proving systems have been kept too for backward compatibility.
* [`r1cs-core`](r1cs/core): Defines core functionalities for rank-1 constraint systems (the circuit synthesizer). 
* [`r1cs-std`](r1cs/gadgets/std): This crate contains elementary "standard" circuits (or, "gadgets"): Boolean operations, native field and elliptic curve arithmetics.  
* [`r1cs-crypto`](r1cs/gadgets/crypto): Provides the circuits for various cryptographic primitives, such as the Poseidon hash, signature schemes, and SNARK verifiers.

In addition, there is a  [`bench-utils`](bench-utils) crate which contains an infrastructure for benchmarking, including macros for timing code segments. 

## Release Note
However, it does not yet serve proof composition. The [`proof-systems`](proof-systems) subcrate [`darlin`](proof-systems/src/darlin) prepares for the full Darlin protocol suite by providing the traits and structs necessary for proof carrying data, and puts simple Coboundary Marlin proofs from [`marlin`](https://github.com/HorizenLabs/marlin) into this framework. It further contains additional tools for scaling verification:
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
More detailed information can be found in our [build guide](BUILD.md).
