#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cargo $CARGOARGS test --all-features
