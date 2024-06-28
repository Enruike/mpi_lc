all: lc_mpi.x
GCC = gcc-9#gcc working version
SRC_DIR := src
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

FILES = conf.c func.c initial_bulk.c initial_channel.c initial_coaxialcyl.c initial_drop.c initial_ellip.c\
		initial_halfcylinder.c initial_halfdrop.c initial_quartercylinder.c initial_quarterdrop.c initial_sandwich.c\
		initial.c main.c output.c relax.c scatter.c energy.c initial_cylinder.c initial_nano_channel.c\
		initial_shell.c initial_not_evolving_shell.c

HEADERS = finite.h

OBJS = read_param.o

WARNS = -diag-disable=10441

lc_mpi.x: $(FILES) $(HEADERS) read_param.o
	mpiicc -gcc-name=$(GCC) -O2 $(FILES) $(OBJS) $(HEADERS) -o lc_mpi.x $(WARNS)
	
debug: $(FILES) $(HEADERS) read_param.o
	mpiicc -gcc-name=$(GCC) -O2 $(FILES) $(OBJS) $(HEADERS) -o lc_mpi.x $(WARNS) -g

read_param.o: read_param.c read_param.h
	icc -gcc-name=$(GCC) -c read_param.c $(WARNS)

move: 
	mv lc_mpi.x ~/Oblates/nanochannel/

test:
	$(OBJ)

clean:
	rm -f $(OBJS) lc_mpi.x
