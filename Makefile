MAIN	= include
CC	= g++
FLAGS	= -Wall -ansi -pedantic -s -O3 -std=c++11 -iquote$(MAIN) -static
DFLAGS	= -Wall -ansi -pedantic -O0 -g -std=c++11 -iquote$(MAIN) -static

SOURCE1	= src/hmm.cpp
BIN	= bin
PROG	= rHMM

.cpp.o:; $(CC) -c $(FLAGS) $<

run:
	mkdir -p $(BIN) && $(CC) -o $(BIN)/$(PROG) $(FLAGS) $(SOURCE1)

debug:
	mkdir -p $(BIN) && $(CC) -o $(BIN)/debug.$(PROG) $(DFLAGS) $(SOURCE1)

clean:
	rm -f $(BIN)/$(PROG)
	rm -f $(BIN)/debug.$(PROG)
