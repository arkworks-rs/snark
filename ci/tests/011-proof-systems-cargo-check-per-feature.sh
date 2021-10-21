#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd proof-systems
cargo $CARGOARGS check --features "gm17" || retval="$?"
cargo $CARGOARGS check --features "groth16" || retval="$?"
cargo $CARGOARGS check --features "darlin" || retval="$?"
exit "$retval"
