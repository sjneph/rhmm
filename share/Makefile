MAIN	= ../include
CC	= g++
FLAGS	= -Wall -ansi -pedantic -s -O3 -iquote$(MAIN) -static -std=c++11
DFLAGS	= -Wall -ansi -pedantic -O0 -g -iquote$(MAIN) -static -std=c++11

SOURCE1	= test1.cpp
TESTBIN	= ../bin
PROG	= test1

.cpp.o:; $(CC) -c $(FLAGS) $<

run:
	mkdir -p $(TESTBIN) && $(CC) -o $(TESTBIN)/$(PROG) $(FLAGS) $(SOURCE1)

debug:
	mkdir -p $(TESTBIN) && $(CC) -o $(TESTBIN)/debug.$(PROG) $(DFLAGS) $(SOURCE1)

clean:
	rm -rf $(TESTBIN)
