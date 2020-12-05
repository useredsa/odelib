# Usage ───────────────────────────────────────────────────────────────────────

#
# 'make'                   compiles all scripts
# 'make clean'             removes all the files generated during compilation
# 'make cleandep'          deletes all .d files
# 'make info'              prints makefile debugging information
#

# Folder Structure ────────────────────────────────────────────────────────────

TARGETSDIR=scripts#        main files folder
INCDIR=include#            headers folder
SRCDIR=src#                source folder/s and file extension
BUILDDIR=build#            object files folder
BINDIR=bin#                binaries folder

# Compilation Process  ────────────────────────────────────────────────────────

CXX=g++-10
INC=-I lib/ -I $(INCDIR)
CXXFLAGS=-std=c++20 -Wall -Werror -Wno-unused -O3 -march=native -mtune=native \
         --fast-math -D EIGEN_DONT_VECTORIZE -D NDEBUG #-ffast-math

# File Sources ────────────────────────────────────────────────────────────────

# patsubst reference:
# https://www.gnu.org/software/make/manual/html_node/Text-Functions.html
# substitution reference: $(var:string1=string2)
# https://www.gnu.org/software/make/manual/html_node/Substitution-Refs.html#Substitution-Refs

# C++ source files that act as scripts (contain a main)
TARGETS=$(shell find $(TARGETSDIR) -type f -name "*.cpp")
# Corresponding binaries
BINS=$(patsubst $(TARGETSDIR)/%,$(BINDIR)/%,$(TARGETS:.cpp=))

# C++ source files for the headers that have them
SRCS=$(shell find $(SRCDIR) -type f -name "*.cpp")
# Corresponding C++ object files 
OBJS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRCS:.cpp=.o))

# Auto-generated dependency files
DEPS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRCS:.cpp=.d)) \
     $(patsubst $(TARGETSDIR)/%,$(BUILDDIR)/%,$(TARGETS:.cpp=.d))

# Rules ───────────────────────────────────────────────────────────────────────

.PHONY: clean cleandep info

all: $(BINS)

clean:
	@rm -rf $(BUILDDIR) $(BINDIR)

cleandep:
	@rm -f $(DEPS)

info:
	@echo "[*] Bin dir:         $(BINDIR)  "
	@echo "[*] Build dir:       ${BUILDDIR}"
	@echo "[*] Targets:         $(TARGETS) "
	@echo "[*] Sources:         $(SRCS)    "
	@echo "[*] Objects:         ${OBJS}    "
	@echo "[*] Dependencies:    ${DEPS}	   "

$(BINDIR)/%: $(TARGETSDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC) $< -o $(@:.d=.o)

# rules to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)

$(BUILDDIR)/%.d: $(TARGETSDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC) $< -MM -MT $(patsubst $(BUILDDIR)/%,$(BINDIR)/%,$(@:.d=)) > $@

$(BUILDDIR)/%.d: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC) $< -MM -MT $(@:.d=.o) > $@

# Include dependencies files

-include $(DEPS)

