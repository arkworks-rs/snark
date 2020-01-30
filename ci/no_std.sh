#!/bin/bash

set -ex

cd algebra
cross build -p algebra --no-default-features --target thumbv6m-none-eabi
cross check --examples -p algebra --no-default-features --target thumbv6m-none-eabi
cd ..

cd r1cs-core
cross build -p r1cs-core --no-default-features --target thumbv6m-none-eabi
cross check --examples -p r1cs-core --no-default-features --target thumbv6m-none-eabi
cd ..
