all: lc_mpi.x

FILES = conf.c func.c initial_bulk.c initial_channel.c initial_coaxialcyl.c initial_drop.c initial_ellip.c initial_halfcylinder.c\
		initial_halfdrop.c initial_quartercylinder.c initial_quarterdrop.c initial_sandwich.c initial.c main.c output.c\
		relax.c scatter.c energy.c initial_cylinder.c

HEADERS = finite.h

OBJS = read_param.o

#WARNS = -diag-disable=10441

lc_mpi.x: $(FILES) $(HEADERS) read_param.o
	mpiicc -O2 $(FILES) $(OBJS) $(HEADERS) -o lc_mpi.x $(WARNS)

read_param.o: read_param.c read_param.h
	icc -c read_param.c

move: 
	mv lc_mpi.x ~/mpi_test/
	
clean:
	rm -f $(OBJS) lc_mpi.x