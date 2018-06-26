SHELL=/bin/sh

SRCS= \
RWeqnSolver.cpp \
Rin.cpp \
Rout.cpp \
Rc.cpp \
system.cpp


HDRS= \
RWeqnSolver.h \
Rin.h \
Rout.h \
Rc.h \
system.h

MAKEFILE=makefile

COMMAND=run

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -O3 -lgsl -lgslcblas

WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas -O3
LIBS=

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

Rin.o : Rin.cpp Rin.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c Rin.cpp -o Rin.o

Rout.o : Rout.cpp Rout.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c Rout.cpp -o Rout.o

Rc.o : Rc.cpp Rc.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c Rc.cpp -o Rc.o

system.o : system.cpp system.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c system.cpp -o system.o

RWeqnSolver.o : RWeqnSolver.cpp RWeqnSolver.h Rin.h Rout.h Rc.h system.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c RWeqnSolver.cpp -o RWeqnSolver.o

clean:
	rm -f $(OBJS) $(COMMAND)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
