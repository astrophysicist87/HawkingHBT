SHELL=/bin/sh

SRCS= \
Zn.cpp \
gauss_quadrature.cpp

HDRS= \
gauss_quadrature.h \
Zn.h \
chebyshev_library.h 

MAKEFILE=makefile

COMMAND=run_Zn

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -pg -g -O3
WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include
 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

#Zn.o : Zn.cpp Zn.h
#	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c Zn.cpp -o Zn.o

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
