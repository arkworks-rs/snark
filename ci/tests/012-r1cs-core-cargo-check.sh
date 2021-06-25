#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cd r1cs/core
cargo $CARGOARGS check
