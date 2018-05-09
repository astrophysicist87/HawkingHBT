SHELL=/bin/sh

SRCS= \
main.cpp \
Zn.cpp \
Yrt.cpp \
gauss_quadrature.cpp

HDRS= \
gauss_quadrature.h \
Zn.h \
Yrt.h \
chebyshev_library.h 

MAKEFILE=makefile

COMMAND=run_Zn_and_Yrt

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -pg -g -O3
WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include
 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
