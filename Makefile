CLINGOROOT=/home/wv/opt/clingo-banane

CXX=clang++-3.8
CXXFLAGS=-std=c++11 -W -Wall -g -O0
CPPFLAGS=-I$(CLINGOROOT)/include
LDFLAGS=-L$(CLINGOROOT)/lib -Wl,-rpath=$(CLINGOROOT)/lib
LDLIBS=-lclingo

TARGET=propagator
SOURCE=main.cpp

OBJECT=$(patsubst %,%.o,$(basename $(SOURCE)))
DEPEND=$(patsubst %,%.d,$(basename $(SOURCE)))

$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $@ $^

%.o: %.cpp
%.o: %.cpp %.d
	$(CXX) -c -MT $@ -MMD -MP -MF $*.Td $(CXXFLAGS) $(CPPFLAGS) $<
	mv -f $*.Td $*.d

.PHONY: clean
format:
	clang-format-3.8 -style="{BasedOnStyle: llvm, IndentWidth: 4, SortIncludes: false, ColumnLimit: 999, AccessModifierOffset: -4}" -i $(SOURCE)

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJECT) $(DEPEND) $(patsubst %,%.Td,$(basename $(SOURCE)))

.PRECIOUS: %.d
%.d: ;

include $(wildcard $(DEPEND))
