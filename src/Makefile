#
# Makefile for nets-nodegroups
#

MAIN = nodegroups
DEPH = group.h
DEPCPP = group.cpp


#
# Global configuration for SNAP makefiles
#

GLIB = glib-core
SNAP = snap-core
GLIBADV = glib-adv
SNAPADV = snap-adv
SNAPEXP = snap-exp

CGLIB = ../snap/$(GLIB)
CSNAP = ../snap/$(SNAP)

EXGLIB = ../snap/$(GLIB)
EXSNAP = ../snap/$(SNAP)
EXGLIBADV = ../snap/$(GLIBADV)
EXSNAPADV = ../snap/$(SNAPADV)
EXSNAPEXP = ../snap/$(SNAPEXP)

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  # Linux flags
  CC = g++
  CXXFLAGS += -std=c++98 -Wall
  CXXFLAGS += -O3 -DNDEBUG
  # turn on for crash debugging, get symbols with <prog> 2>&1 | c++filt
  #CXXFLAGS += -g -rdynamic
  CXXOPENMP = -fopenmp
  LDFLAGS +=
  LIBS += -lrt

else ifeq ($(UNAME), Darwin)
  # OS X flags
  CC = g++
  CXXFLAGS += -std=c++98 -Wall
  CXXFLAGS += -O3 -DNDEBUG
  CXXOPENMP = -fopenmp
  CLANG := $(shell g++ -v 2>&1 | grep clang | cut -d " " -f 2)
  ifeq ($(CLANG), LLVM)
    CXXOPENMP =
    CXXFLAGS += -DNOMP
  endif

  LDFLAGS +=
  LIBS +=

else ifeq ($(shell uname -o), Cygwin)
  # Cygwin flags
  CC = g++
  CXXFLAGS += -Wall -D__STDC_LIMIT_MACROS
  CXXFLAGS += -O3 -DNDEBUG
  CXXOPENMP = -fopenmp
  LDFLAGS +=
  LIBS +=

endif


#
# Makefile for non-Microsoft compilers
#

all: $(MAIN)
	make -C experiments

# COMPILE
$(MAIN): $(MAIN).cpp $(DEPH) $(DEPCPP) $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(MAIN) $(MAIN).cpp $(DEPCPP) $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(EXSNAP)/Snap.o: 
	make -C $(EXSNAP)

clean:
	make clean -C experiments
	rm -f *.o $(MAIN) $(MAIN).exe *.user
	rm -rf Debug Release
