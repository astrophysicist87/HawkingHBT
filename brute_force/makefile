SHELL=/bin/sh

SRCS= \
RWeqnSolver.cpp

HDRS=

MAKEFILE=makefile

COMMAND=run

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -O3 -lgsl -lgslcblas

WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas -O3
LIBS=

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

RWeqnSolver.o : RWeqnSolver.cpp
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c RWeqnSolver.cpp -o RWeqnSolver.o

clean:
	rm -f $(OBJS) $(COMMAND)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
