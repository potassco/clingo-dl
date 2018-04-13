undefine CXX
-include FLAGS

CLINGOROOT?=/home/wv/opt/clingo-banane
LUA_INCLUDE_DIR?=/usr/include/lua5.3

CXX?=c++
CXXFLAGS?=-std=c++14 -W -Wall -O3 -DNDEBUG
CPPFLAGS?=-I$(CLINGOROOT)/include
LDFLAGS?=-L$(CLINGOROOT)/lib -Wl,-rpath=$(CLINGOROOT)/lib
LDLIBS?=-lclingo

TARGET=clingoDL
LIBTARGET=lib$(TARGET)
LIBLUATARGET=liblua$(TARGET)
SOURCE=main.cpp
OBJECT=$(patsubst %,%.o,$(basename *))
DEPEND=$(patsubst %,%.d,$(basename $(SOURCE)))

all: $(TARGET) $(LIBTARGET) $(LIBLUATARGET)

$(TARGET): main.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)
main.o: main.cpp
	$(CXX) -c $^ -MT $@ -MMD -MP -MF $(CXXFLAGS) $(CPPFLAGS)

$(LIBTARGET): libmain.o
	$(CXX) -o $@.so $^ $(LDFLAGS) -shared $(LDLIBS)
libmain.o: main.cpp
	$(CXX) -fPIC -c $^ -o $@ -MT $@ -MMD -MP -MF $(CXXFLAGS) $(CPPFLAGS) 

$(LIBLUATARGET): libluamain.o
	$(CXX) -o $@.so $^ $(LDFLAGS) -shared $(LDLIBS)
libluamain.o: main.cpp
	$(CXX) -fPIC -c $^ -o $@ -MT $@ -MMD -MP -MF $(CXXFLAGS) $(CPPFLAGS) -DWITH_LUA -I$(LUA_INCLUDE_DIR)

FLAGS:
	rm -f FLAGS
	echo "CXX:=$(CXX)" >> FLAGS
	echo "CXXFLAGS:=$(CXXFLAGS)" >> FLAGS
	echo "CPPFLAGS:=$(CPPFLAGS)" >> FLAGS
	echo "LDFLAGS:=$(LDFLAGS)" >> FLAGS
	echo "LDLIBS:=$(LDLIBS)" >> FLAGS
	echo "LUA_INCLUDE_DIR:=$(LUA_INCLUDE_DIR)" >> FLAGS

.PHONY: format
format:
	clang-format-3.8 -style="{BasedOnStyle: llvm, IndentWidth: 4, SortIncludes: false, ColumnLimit: 256, AccessModifierOffset: -4, BreakBeforeBraces: Custom, BraceWrapping: {BeforeElse: true, BeforeCatch: true}, BreakConstructorInitializersBeforeComma: true, AlwaysBreakTemplateDeclarations: true, AlignAfterOpenBracket: AlwaysBreak, AllowShortBlocksOnASingleLine: true, IndentCaseLabels: true}" -i $(SOURCE)

.PHONY: clean
clean:
	rm -f $(TARGET) $(LIBTARGET).so $(LIBLUATARGET).so $(OBJECT) $(DEPEND) $(patsubst %,%.Td,$(basename $(SOURCE)))

.PRECIOUS: %.d
%.d: ;

include $(wildcard $(DEPEND))
