## Pending
- [\#30](https://github.com/arkworks-rs/groth16/pull/30) Add proof input preprocessing.

### Breaking changes
- [\#21](https://github.com/arkworks-rs/groth16/pull/21) Change the `generate_parameters` interface to take generators as input.

### Features

### Improvements

### Bug fixes

## v0.2.0

### Breaking changes
- [\#4](https://github.com/arkworks-rs/groth16/pull/4) Change groth16's logic to implement the `SNARK` trait.
- Minimum version on crates from `arkworks-rs/algebra` and `arkworks-rs/curves` is now `v0.2.0`
- [\#24](https://github.com/arkworks-rs/groth16/pull/24) Switch from `bench-utils` to `ark_std::perf_trace`

### Features
- [\#5](https://github.com/arkworks-rs/groth16/pull/5) Add R1CS constraints for the groth16 verifier.
- [\#8](https://github.com/arkworks-rs/groth16/pull/8) Add benchmarks for the prover
- [\#16](https://github.com/arkworks-rs/groth16/pull/16) Add proof re-randomization

### Improvements
- [\#9](https://github.com/arkworks-rs/groth16/pull/9) Improve memory consumption by manually dropping large vectors once they're no longer needed

### Bug fixes
- [c9bc5519](https://github.com/arkworks-rs/groth16/commit/885b9b569522f59a7eb428d1095f442ec9bc5519) Fix parallel feature flag
- [\#22](https://github.com/arkworks-rs/groth16/pull/22) Compile with `panic='abort'` in release mode, for safety of the library across FFI boundaries.

## v0.1.0

_Initial release_
