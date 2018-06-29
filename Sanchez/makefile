SHELL=/bin/sh

SRCS= \
SanchezAsymptotics.cpp

HDRS=

MAKEFILE=makefile

COMMAND=run_Sanchez

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -g -O3
WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include
 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
