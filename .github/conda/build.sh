#!/bin/bash

set -euo pipefail

cmake -G Ninja -S . -B build \
    -DCLINGODL_MANAGE_RPATH=Off \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DCMAKE_INSTALL_LIBDIR="lib" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DPYCLINGODL_ENABLE="require" \
    -DPython_EXECUTABLE="${PYTHON}"

cmake --build build --verbose
cmake --build build --target install --verbose
