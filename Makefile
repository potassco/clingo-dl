BUILD_TYPE:=debug
POTASSCO_PREFIX:=${HOME}/.local/opt/potassco/$(BUILD_TYPE)
CXXFLAGS=-Wall -Wextra -Wpedantic -Werror
define cmake_options
-S . \
-B "build/$(BUILD_TYPE)" \
-DCMAKE_INSTALL_PREFIX="$(POTASSCO_PREFIX)" \
-DCMAKE_CXX_FLAGS="$(CXXFLAGS)" \
-DCLINGODL_BUILD_TESTS=On \
-DCMAKE_EXPORT_COMPILE_COMMANDS=On
endef

ifeq ($(BUILD_TYPE),profile)
	cmake_options += -DCMAKE_BUILD_TYPE="RelWithDebInfo" -DCLINGODL_PROFILE=On
else
	cmake_options += -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)"
endif


.PHONY: all test compdb configure

all: configure
	$(MAKE) -C "build/$(BUILD_TYPE)"

test: all
	$(MAKE) -C "build/$(BUILD_TYPE)" test

%: configure
	@TERM=dumb MAKEFLAGS= MFLAGS= cmake --build "build/$(BUILD_TYPE)" --target "$@"

compdb:
	compdb -p "build/$(BUILD_TYPE)" list -1 > compile_commands.json

configure: build/$(BUILD_TYPE)/Makefile

build/$(BUILD_TYPE)/Makefile:
	cmake $(cmake_options)

Makefile:
	:
