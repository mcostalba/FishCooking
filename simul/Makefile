EXE = simul

OBJS = simul.o

CXXFLAGS = -g -Wall -Wextra -Wshadow -Wcast-qual -ansi -pedantic -std=c++11 #-fno-rtti -03

$(EXE): $(OBJS)
	$(CXX) -o $@ $(OBJS) --std=c11

clean:
	$(RM) $(EXE) *.o .depend test

.depend:
	$(CXX) -MM -std=c++11 $(OBJS:.o=.cc) >$@

-include .depend
