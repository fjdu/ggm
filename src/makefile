cpl?=gfortran

cflags=-c
exe_name?=mo

all: lflags=-O3 -fbounds-check -o
all: $(exe_name)

objects = common_MC_MO.o common_MC_MO_TP.o subroutines_share_triv.o mo.o subroutines_share_proc.o subroutines_mo.o opkd*.o

$(exe_name): $(objects)
	$(cpl) $(lflags) $(exe_name) $(objects)

mo.o: mo.f90 common_MC_MO.o common_MC_MO_TP.o
	$(cpl) $(cflags) mo.f90

common_MC_MO.o: common_MC_MO.f90
	$(cpl) $(cflags) common_MC_MO.f90

common_MC_MO_TP.o: common_MC_MO_TP.f90
	$(cpl) $(cflags) common_MC_MO_TP.f90

subroutines_share_triv.o: subroutines_share_triv.f90
	$(cpl) $(cflags) subroutines_share_triv.f90

subroutines_share_proc.o: subroutines_share_proc.f90
	$(cpl) $(cflags) subroutines_share_proc.f90

subroutines_mo.o: subroutines_mo.f90
	$(cpl) $(cflags) subroutines_mo.f90

opkd*.o: opkd*.f
	$(cpl) $(cflags) opkd*.f

clean:
	rm *.o *.mod

to193:
	rsync -acvz ~/_c/DeuterationGasGrain 193:~/code/
