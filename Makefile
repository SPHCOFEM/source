SRC = errutils.c ioutils.c allocator.c kernel.c volume.c nns.c move.c interpolation.c accelerations.c rigid_bodies.c matrices.c spatial_rotation.c contact_forces.c contacts.c matrices_operations.c fem.c fem_bar.c fem_beam.c fem_triangle.c fem_quad.c fem_membrane.c fem_shell.c fem_tetrahedron.c cgm.c store.c freeall.c energy.c time_step.c user_defined_material.c sphcofem.c
OBJ = errutils.o ioutils.o allocator.o kernel.o volume.o nns.o move.o interpolation.o accelerations.o rigid_bodies.o matrices.o spatial_rotation.o contact_forces.o contacts.o matrices_operations.o fem.o fem_bar.o fem_beam.o fem_triangle.o fem_quad.o fem_membrane.o fem_shell.o fem_tetrahedron.o cgm.o store.o freeall.o energy.o time_step.o user_defined_material.o sphcofem.o

#CFLG = -g -Wall -O2 -D__SDIR__='"./"' -fno-math-errno
CFLG = -g -Wall -O2 -D__SDIR__='"./"'
#CFLG = -g3 -O2
CC = gcc
#CC = cc

TARGET = sphcofem

.c.o:
	$(CC) -c $(CFLG) $<

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) -lm
