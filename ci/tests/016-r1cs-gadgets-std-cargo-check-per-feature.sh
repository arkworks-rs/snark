#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd r1cs/gadgets/std
cargo $CARGOARGS check --features "full" || retval="$?"
cargo $CARGOARGS check --features "bls12_377" || retval="$?"
cargo $CARGOARGS check --features "bn_382" || retval="$?"
cargo $CARGOARGS check --features "edwards_bls12" || retval="$?"
cargo $CARGOARGS check --features "edwards_sw6" || retval="$?"
cargo $CARGOARGS check --features "jubjub" || retval="$?"
cargo $CARGOARGS check --features "mnt4_753" || retval="$?"
cargo $CARGOARGS check --features "mnt6_753" || retval="$?"
cargo $CARGOARGS check --features "tweedle" || retval="$?"
cargo $CARGOARGS check --features "nonnative" || retval="$?"
cargo $CARGOARGS check --features "nonnative, density-optimized" || retval="$?"
exit "$retval"
