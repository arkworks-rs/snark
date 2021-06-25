#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cd proof-systems
cargo $CARGOARGS check --all-features --benches --examples
