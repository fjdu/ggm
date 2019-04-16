CPL?=g++
CPLOPT=-std=c++11
FCC?=gfortran
LINKER?=gfortran
LINKOPT=-lc++
#  -pedantic -Wall -Wextra

target?=re

all: $(target)

OBJS=main.o rate_equation_lsode.o ode_wrapper.o types.o constants.o opkdmain.o opkda1.o opkda2.o logistics.o calculate_reaction_rate.o utils.o

$(target): $(OBJS)
	$(LINKER) $(LINKOPT) -o $(target) $(OBJS)

types.o: types.hpp types.cpp
	$(CPL) $(CPLOPT) -c types.cpp

constants.o: constants.cpp
	$(CPL) $(CPLOPT) -c constants.cpp

main.o: main.cpp rate_equation_lsode.hpp logistics.hpp types.hpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c main.cpp

logistics.o: logistics.cpp
	$(CPL) $(CPLOPT) -c logistics.cpp

calculate_reaction_rate.o: calculate_reaction_rate.cpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c calculate_reaction_rate.cpp

rate_equation_lsode.o: rate_equation_lsode.cpp
	$(CPL) $(CPLOPT) -c rate_equation_lsode.cpp

utils.o: utils.cpp
	$(CPL) $(CPLOPT) -c utils.cpp

ode_wrapper.o: ode_wrapper.f90
	$(FCC) -c ode_wrapper.f90

opkdmain.o: opkdmain.f opkda1.f opkda2.f
	$(FCC) -c opkd*.f

clean:
	rm *.o *.gch
