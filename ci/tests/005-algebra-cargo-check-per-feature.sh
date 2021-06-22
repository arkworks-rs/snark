#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd algebra
cargo $CARGOARGS check --features "parallel" || retval="$?"
cargo $CARGOARGS check --features "fft" || retval="$?"
cargo $CARGOARGS check --features "n_fold" || retval="$?"
cargo $CARGOARGS check --features "bls12_377" || retval="$?"
cargo $CARGOARGS check --features "bls12_381" || retval="$?"
cargo $CARGOARGS check --features "edwards_bls12" || retval="$?"
cargo $CARGOARGS check --features "edwards_sw6" || retval="$?"
cargo $CARGOARGS check --features "jubjub" || retval="$?"
cargo $CARGOARGS check --features "sw6" || retval="$?"
cargo $CARGOARGS check --features "mnt4_753" || retval="$?"
cargo $CARGOARGS check --features "mnt6_298" || retval="$?"
cargo $CARGOARGS check --features "mnt6_753" || retval="$?"
cargo $CARGOARGS check --features "bn_382" || retval="$?"
cargo $CARGOARGS check --features "tweedle" || retval="$?"
cargo $CARGOARGS check --features "full" || retval="$?"
exit "$retval"
