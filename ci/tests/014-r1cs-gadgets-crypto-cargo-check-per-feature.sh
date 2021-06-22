#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd r1cs/gadgets/crypto
cargo $CARGOARGS check --features "commitment" || retval="$?"
cargo $CARGOARGS check --features "merkle_tree" || retval="$?"
cargo $CARGOARGS check --features "prf" || retval="$?"
cargo $CARGOARGS check --features "signature" || retval="$?"
cargo $CARGOARGS check --features "vrf" || retval="$?"
cargo $CARGOARGS check --features "nizk" || retval="$?"
cargo $CARGOARGS check --features "mnt4_753" || retval="$?"
cargo $CARGOARGS check --features "mnt6_753" || retval="$?"
cargo $CARGOARGS check --features "bn_382" || retval="$?"
cargo $CARGOARGS check --features "tweedle" || retval="$?"
exit "$retval"
