# ----------------------------------------------------------------
# environment
CC		= gcc
FC		= 

# ----------------------------------------------------------------
# options

CFLAGS		= -O3 -fopenmp
FFLAGS		= 

# ----------------------------------------------------------------
# sources and objects

C_SRC		= cpi_m.c
F_SRC		= 

C_OBJ		= $(C_SRC:.c=)
F_OBJ		= $(F_SRC:.f=)

# ----------------------------------------------------------------
# executables

EXEC		= $(C_OBJ) 

all:		$(EXEC)

$(C_OBJ):	$(C_SRC)
	$(CC) -o $@ $(CFLAGS) $(C_SRC) -lm


# ----------------------------------------------------------------
# rules

.c.:
	$(CC) -o $* $(CFLAGS) -c $<

.f.:
	$(FC) -o $* $(FFLAGS) -c $<

# ----------------------------------------------------------------
# clean up

clean:
	/bin/rm -f $(EXEC) $(F_SRC:.f=.o)

# ----------------------------------------------------------------
# End of Makefile
