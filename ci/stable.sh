#!/bin/bash

set -ex

cargo test --all -- --skip dpc --skip integration_test
cargo check --examples --all
