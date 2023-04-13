all: lc_mpi.x

FILES = conf.c func.c initial_bulk.c initial_channel.c initial_coaxialcyl.c initial_drop.c initial_ellip.c initial_halfcylinder.c\
		initial_halfdrop.c initial_quartercylinder.c initial_quarterdrop.c initial_sandwich.c initial.c main.c output.c\
		read_param.c relax.c scatter.c energy.c initial_cylinder.c

HEADERS = finite.h

lc_mpi.x: $(FILES) $(HEADERS)
	mpiicc -O2 $(FILES) -o lc_mpi.x