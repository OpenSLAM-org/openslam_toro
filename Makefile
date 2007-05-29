all: doall

################################

TARGETS =  libtoro.a toro

libtoro.a: posegraph2.o treeoptimizer2.o

toro: posegraph2.o treeoptimizer2.o  toro.o

################################
#
# CONFIG
#

SILENT = @
ECHO = @echo
MAKE = make

CXX=/usr/bin/g++
LINK = $(CXX)
AR = ar cr
RM = rm -rf
RANLIB = ranlib

ifdef CYGWIN
CXXFLAGS += -DCYGWIN
endif

CXXFLAGS += -O3 -ffast-math -Wall -W 
IFLAGS =
LFLAGS =


################################
#
# RULES
#


banner = $(SILENT) \
	echo ; \
	echo "*******************************************************************" ; \
	echo "* TORO (c) by Giorgio Grisetti, Cyrill Stachniss, Wolfram Burgard *" ; \
	echo "*******************************************************************" ; \


doall:
	$(banner)
	$(SILENT) for i in $(filter %.a, $(TARGETS)) xxxx ; do \
		if test ! "$$i" = "xxxx" ; then \
			if ! $(MAKE) $$i ; then \
                                exit -1; \
                        fi; \
		fi ; \
	done 
	$(SILENT) for i in $(filter-out %.a, $(TARGETS)) xxxx ; do \
		if test ! "$$i" = "xxxx" ; then \
			if ! $(MAKE) $$i ; then \
                                exit -1; \
                        fi; \
		fi ; \
	done 

clean:
	$(banner)
	$(ECHO) " ---- Cleaning up "
	$(SILENT) $(RM) *.o *.a *.so *.exe core a.out Makefile.depend Makefile.depend.bak $(TARGETS) 

%.o:	%.cpp
	$(ECHO) " ---- Compiling $< to $@ (CXX)"
	$(SILENT) $(CXX) $(CXXFLAGS) $(IFLAGS) -c $< 

%.a:	
	$(ECHO) " ---- Archiving $^ into $@ (CXX)"
	$(SILENT) $(AR) $@ $^
	$(SILENT) $(RANLIB) $@

%:	%.o
	$(ECHO) " ---- Linking $^ to $@ (CXX)"
	$(SILENT) $(LINK) $(CXXFLAGS) $(IFLAGS) $(filter %.o, $^) $(filter %.a, $^) -o $@ -L. $(patsubst lib%.a,-l%,$(filter %.a, $^)) $(LFLAGS) $(LFLAGS_POST)
