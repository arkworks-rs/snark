<h1 align="center">ark-abpr21 (RO_Based)</h1>

<p align="center">
    <img src="https://github.com/arkworks-rs/snark/workflows/CI/badge.svg?branch=master">
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/groth16/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/arkworks-rs/snark"><img src="https://deps.rs/repo/github/arkworks-rs/snark/status.svg"></a>
</p>

The arkworks ecosystem consist of Rust libraries for designing and working with __zero knowledge succinct non-interactive arguments (zkSNARKs)__. This repository contains an efficient Rust implementation of the RO-based simulation extractable variant of [[Groth16]](https://eprint.iacr.org/2016/260) zkSNARK presented in Section 4 of [[ABPR21]](https://eprint.iacr.org/2020/1306), implemented by Oussama Amine (University of Oslo) and Karim Baghery (KU Leuven).

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

Below we compare the performance of several weakly- and strongly-simulation-extractable zkSNARKs. Note that Groth's zk-SNARK is proven to achieve weak simulation-extractability [[BKSV20]](https://eprint.iacr.org/2020/811). We benchmark the zkSNARKs on an R1CS instance for different curves and report proving and verifying times for each. 
All experiments were performed on a desktop machine running Ubuntu 20.4.2 LTS and equipped with an Intel Core i9-9900 processor at base frequency 3.1 GHz, and 128GB of memory.
Proving was performed in multi-threaded mode, with 16 threads, while proof verification was tested in single-threaded mode.

Key: _SE_ = Simulation-extractable, _RO_ = Random Oracle, _CRH_ = Collision Resistant Hash.

| curve | zkSNARK | security | per-constraint proving time (ns) | proof size (bytes) |  verifying time, 1 proof (ms) |  verifying time, 100 proofs (s) | verifying time, 1000 proofs (s) | 
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | 
| BLS12-381 | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 5026  | 127.5 | 1.90 |  0.19  |  1.90|
| BLS12-381 | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 11042 | 127.5 | 3.32 |  0.322 |  3.32|
| BLS12-381 | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 5052  | 223.1 | 3.52 |  0.352 |  3.52|
| BLS12-381 | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)   | Strong SE | 5042  | 223.1 | 4.85 |  0.360 |  3.50|
| BLS12-381 | ABPR21-RO                                                        | Strong SE | 5041  | 191.2 | 2.39 |  0.194 |  1.91|
| MNT4-298  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 4830  | 149.0 | 2.67 |  0.267 |  2.67|
| MNT4-298  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 10025 | 149.0 | 3.80 |  0.380 |  3.80|
| MNT4-298  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 4879  | 260.7 | 4.32 |  0.432 |  4.32|
| MNT4-298  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)   | Strong SE | 4881  | 260.7 | 4.45 |  0.311 |  3.05|
| MNT4-298  |  ABPR21-RO                                                       | Strong SE | 4875  | 223.5 | 3.33 |  0.271 |  2.68|
| MTN6-298  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 5794  | 186.2 | 4.94 |  0.494  |  4.91|
| MTN6-298  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 11427 | 186.2 | 7.07 |  0.707 |  7.07|
| MTN6-298  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 5831  | 335.2 | 8.07 |  0.807 |  8.07|
| MTN6-298  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)   | Strong SE | 5824  | 335.2 | 8.34 |  0.582 |  5.72|
| MTN6-298  | ABPR21-RO                                                        | Strong SE | 5810  | 298.0 | 6.11 |  0.501 |  4.97|
| MNT4-753  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 30247 | 376.5 | 29.1 |  2.91 |  29.1|
| MNT4-753  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 83120 | 376.5 | 41.6 |  4.16 |  41.6|
| MNT4-753  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 30863 | 658.8 | 47.3 |  4.73 |  47.3|
| MNT4-753  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)   | Strong SE | 30887 | 658.8 | 45.5 |  3.41 |  33.8|
| MNT4-753  | ABPR21-RO                                                        | Strong SE | 30760 | 564.7 | 33.9 |  2.94 |  29.2|
| MTN6-753  | [Gro16](https://github.com/arkworks-rs/groth16)                  | Weak SE   | 33298  | 470.6 | 53.6 |  5.36  | 53.6|
| MTN6-753  | [GM17](https://github.com/arkworks-rs/gm17)                      | Strong SE | 83121  | 470.6 | 76.9 |  7.69 |  76.9|
| MTN6-753  | [BG18](https://github.com/Baghery/ABPR21/tree/BG18)              | Strong SE | 33358  | 847.1 | 88.5 |  8.85 |  88.5|
| MTN6-753  | [ABPR21-CRH](https://github.com/Baghery/ABPR21/tree/CRH_Based)   | Strong SE | 33359  | 847.1 | 85.4 |  6.33 |  63.1|
| MTN6-753  | ABPR21-RO                                                        | Strong SE | 33345  | 753.0 | 64.4 |  5.42 |  53.8|

## License

This library is licensed under either of the following licenses, at your discretion.

* Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in this library by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## References

[\[ABPR21\]](https://eprint.iacr.org/2020/1306): Simulation Extractable Versions of Groth’s zk-SNARK Revisited. Oussama Amine, Karim Baghery, Zaira Pindado, and Carla Ràfols. 2021. Note: Full version of \[BPR20\].

[\[BPR20\]](https://eprint.iacr.org/2020/1306): Simulation Extractable Versions of Groth’s zk-SNARK Revisited. Karim Baghery, Zaira Pindado, and Carla Ràfols. CANS'20. 2020.