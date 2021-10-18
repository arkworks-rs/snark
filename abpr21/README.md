<h1 align="center">ark-abpr21 (RO_Based)</h1>

<p align="center">
    <img src="https://github.com/arkworks-rs/groth16/workflows/CI/badge.svg?branch=master">
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/arkworks-rs/groth16"><img src="https://deps.rs/repo/github/arkworks-rs/groth16/status.svg"></a>
</p>

The arkworks ecosystem consist of Rust libraries for designing and working with __zero knowledge succinct non-interactive arguments (zkSNARKs)__. This repository contains an efficient Rust implementation of the RO-based simulation extractable variant of [[Groth16]](https://eprint.iacr.org/2016/260) zk-SNARK presented in Section 4 of [[ABPR21]](https://eprint.iacr.org/2020/1306) which is the extended version of the paper [[BPR20]](https://link.springer.com/chapter/10.1007/978-3-030-65411-5_22) appeared in the proceedings of CANS 2020. The imlementations are done by Oussama Amine (University of Oslo) and Karim Baghery (KU Leuven).

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone -b RO_Based https://github.com/Baghery/ABPR21.git
cd RO_Based
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test
```
and for benchmarking the scheme with `RAYON_NUM_THREADS=4` threads, run the following command,  
```bash
RAYON_NUM_THREADS=4 cargo bench --no-default-features --features "std parallel" -- --nocapture
```

## Empirical performance

Below is the empricial performance of several weak and strong simulation extractable zk-SNARKs in `Arkworks`. Note that Groth's zk-SNARK is proven to achieve weak simulation extractability [[BKSV20]](https://eprint.iacr.org/2020/811).  
We benchmark the zk-SNARKs on an R1CS instance for different curves and report proving and verifying times for each constrain with 100 iterations for the prover and 10.000 iterations for the verification. 
All experiments are done on a desktop machine with Ubuntu 20.4.2 LTS, an Intel Core i9-9900 processor at base frequency 3.1 GHz, and 128GB of memory. 
Proof generations are done in the multi-thread mode, with 16 threads, while proof verifications are done in a single-thread mode. In the verification of our constructions, in the case of verifying more than one proof, we use Multi-Scalar Multiplication (MSM) techniques to optimize the computation of exponentiations in $G_2$ and $G_T$. 

Abbreviations used: <i>SE</i> = Simulation Extractable, <i>PCPT</i> = Per-Constraint Proving Time, <i>ns</i> = nanoseconds, <i>RO</i> = Random Oracle, <i>CRH</i> = Collision Resistant Hash.

| Curve | zk-SNARK | Secuiry | PCPT, ns | Proof, bytes |  Verifier, 1 proof |  Verifier, 100 proofs | Verifier, 1000 proofs | 
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | 
| BLS12-381 | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 5026  | 127.5 | 1.90 ms |  0.19 sec   |  1.90 sec |
| BLS12-381 | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 11042 | 127.5 | 3.32 ms |  0.322 sec  |  3.32 sec |
| BLS12-381 | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 5052  | 223.1 | 3.52 ms |  0.352 sec  |  3.52 sec |
| BLS12-381 | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)    | Strong SE | 5042  | 223.1 | 4.85 ms |  0.360 sec  |  3.50 sec |
| BLS12-381 | [ABPR21-RO](https://github.com/Baghery/ABPR21/tree/RO_Based)      | Strong SE | 5041  | 191.2 | 2.39 ms |  0.194 sec  |  1.91 sec |
| **MNT4-298** | [Gro16](https://github.com/arkworks-rs/groth16)               | Weak SE   | 4830  | 149.0 | 2.67 ms |  0.267 sec  |  2.67 sec |
| **MNT4-298** | [GM17](https://github.com/arkworks-rs/gm17)                   | Strong SE | 10025 | 149.0 | 3.80 ms |  0.380 sec  |  3.80 sec |
| **MNT4-298** | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)           | Strong SE | 4879  | 260.7 | 4.32 ms |  0.432 sec  |  4.32 sec |
| **MNT4-298** | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based) | Strong SE | 4881  | 260.7 | 4.45 ms |  0.311 sec  |  3.05 sec |
| **MNT4-298** | [ABPR21-RO](https://github.com/Baghery/ABPR21/tree/RO_Based)   | Strong SE | 4875  | 223.5 | 3.33 ms |  0.271 sec  |  2.68 sec |
| MTN6-298  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 5794  | 186.2 | 4.94 ms |  0.494 sec   |  4.91 sec |
| MTN6-298  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 11427 | 186.2 | 7.07 ms |  0.707 sec  |  7.07 sec |
| MTN6-298  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 5831  | 335.2 | 8.07 ms |  0.807 sec  |  8.07 sec |
| MTN6-298  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)    | Strong SE | 5824  | 335.2 | 8.34 ms |  0.582 sec  |  5.72 sec |
| MTN6-298  | [ABPR21-RO](https://github.com/Baghery/ABPR21/tree/RO_Based)      | Strong SE | 5810  | 298.0 | 6.11 ms |  0.501 sec  |  4.97 sec |
| **MNT4-753** | [Gro16](https://github.com/arkworks-rs/groth16)               | Weak SE   | 30247 | 376.5 | 29.1 ms |  2.91 sec  |  29.1 sec |
| **MNT4-753** | [GM17](https://github.com/arkworks-rs/gm17)                   | Strong SE | 83120 | 376.5 | 41.6 ms |  4.16 sec  |  41.6 sec |
| **MNT4-753** | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)           | Strong SE | 30863 | 658.8 | 47.3 ms |  4.73 sec  |  47.3 sec |
| **MNT4-753** | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based) | Strong SE | 30887 | 658.8 | 45.5 ms |  3.41 sec  |  33.8 sec |
| **MNT4-753** | [ABPR21-RO](https://github.com/Baghery/ABPR21/tree/RO_Based)   | Strong SE | 30760 | 564.7 | 33.9 ms |  2.94 sec  |  29.2 sec |
| MTN6-753  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 33298  | 470.6 | 53.6 ms |  5.36 sec   | 53.6 sec |
| MTN6-753  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 83121  | 470.6 | 76.9 ms |  7.69 sec  |  76.9 sec |
| MTN6-753  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 33358  | 847.1 | 88.5 ms |  8.85 sec  |  88.5 sec |
| MTN6-753  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)    | Strong SE | 33359  | 847.1 | 85.4 ms |  6.33 sec  |  63.1 sec |
| MTN6-753  | [ABPR21-RO](https://github.com/Baghery/ABPR21/tree/RO_Based)      | Strong SE | 33345  | 753.0 | 64.4 ms |  5.42 sec  |  53.8 sec |

## License

This library is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in this library by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## Acknowledgements

This work was supported by:
the Defense Advanced Research Projects Agency (DARPA) under Contract No. HR001120C0085; 
a Google Faculty Award;
the National Science Foundation;
the UC Berkeley Center for Long-Term Cybersecurity;
and donations from the Ethereum Foundation, the Interchain Foundation, and Qtum.

An earlier version of this library was developed as part of the paper *"[ZEXE: Enabling Decentralized Private Computation][zexe]"*.

[zexe]: https://ia.cr/2018/962

