BUILD_TYPE:=debug
CLINGO_DIR:=${HOME}/.local/opt/potassco/$(BUILD_TYPE)/lib/cmake/Clingo
CLANG_TIDY_WARNINGS:=clang-analyzer-*,readability-*,modernize-*,cppcoreguidelines-*,performance-*,bugprone-*,-modernize-use-trailing-return-type,-readability-magic-numbers,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-cppcoreguidelines-pro-type-union-access,-cppcoreguidelines-pro-type-vararg
CLANG_TIDY_ERRORS:=
CLANG_TIDY:=clang-tidy;--extra-arg-before=-stdlib=libc++;-header-filter=.*hh;-checks=$(CLANG_TIDY_WARNINGS);-warnings-as-errors=$(CLANG_TIDY_ERRORS)
CXXFLAGS=-Wall -Wextra -Wpedantic -Werror

.PHONY: all configure compdb

all: configure
	@TERM=dumb MAKEFLAGS= MFLAGS= cmake --build "build/$(BUILD_TYPE)" --target all

test: all
	@cmake --build "build/$(BUILD_TYPE)" --target "test"


%: configure
	@cmake --build "build/$(BUILD_TYPE)" --target "$@"

# compdb can be installed with pip
compdb: configure
	compdb -p "build/$(BUILD_TYPE)" list -1 > compile_commands.json

configure: build/$(BUILD_TYPE)/build.ninja

build/$(BUILD_TYPE)/build.ninja:
	cmake -G Ninja -H. -B"build/$(BUILD_TYPE)" -DCMAKE_INSTALL_PREFIX=${HOME}/.local/opt/potassco/$(BUILD_TYPE) -DCMAKE_CXX_FLAGS="$(CXXFLAGS)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DCMAKE_CXX_CLANG_TIDY:STRING="$(CLANG_TIDY)" -DClingo_DIR="$(CLINGO_DIR)" -DCLINGODL_BUILD_TESTS=On -DCMAKE_EXPORT_COMPILE_COMMANDS=On

Makefile:
	:
