#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cd primitives
cargo $CARGOARGS check --all-features --benches
