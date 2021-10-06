#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd primitives
cargo $CARGOARGS check || retval="$?"
cargo $CARGOARGS check --all-features --tests || retval="$?"
exit "$retval"
