mkdir build

cmake -G "Ninja" -H. -Bbuild ^
    -DCLINGODL_MANAGE_RPATH=Off ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_CXX_COMPILER="%CXX%" ^
    -DCMAKE_INSTALL_BINDIR="." ^
    -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
    -DPYCLINGODL_ENABLE="require" ^
    -DPython_ROOT_DIR="%PREFIX%"

cmake --build build
cmake --build build --target install
