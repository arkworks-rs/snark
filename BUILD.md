# Build guide

The library compiles on the `1.51.0 stable` toolchain of the Rust compiler. To install Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the appropriate Rust toolchain by invoking:
```bash
rustup install 1.51.0
```
After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/HorizenOfficial/ginger-lib.git
cd ginger-lib
cargo build --release
```
This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test --all-features
```
This library comes with several benchmarks.
Some of the benchmarks in [`algebra`](algebra/benches) crate require the nightly Rust toolchain (we suggest to use `nightly-2021-04-25)`; to install this, run `rustup install nightly-2021-04-25`. Then, to run benchmarks, run the following command: 
```bash
cargo +nightly-2021-04-25 bench --all-features 
```
Other benchmarks using [Criterion]() crate are present in [`algebra`](algebra/benches), [`primitives`](primitives/benches/crypto_primitives) and [`proof-systems`](proof-systems) crates and can be ran with `stable` Rust.

Compiling with `adcxq`, `adoxq` and `mulxq` instructions can lead to a 30-70% speedup. These are available on most `x86_64` platforms (Broadwell onwards for Intel and Ryzen onwards for AMD). Run the following command:
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo test/build/bench --features llvm_asm
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
