CPL?=g++
CFLAGS?=-std=c++11 -I/Users/fjdu/Dropbox/codes_from_others/sundials-4.1.0/builddir/include/ -L/Users/fjdu/Dropbox/codes_from_others/sundials-4.1.0/builddir/lib/ -lsundials_cvodes -lsundials_sunlinsolspgmr
#  -pedantic -Wall -Wextra

target?=re

all: $(target)
	$(CPL) $(CFLAGS) main.cpp -o $(target)

$(target): main.cpp

main.cpp: logistics.cpp types.hpp update.hpp
