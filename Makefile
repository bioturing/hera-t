CC = gcc

AR = ar

LIBS = -pthread -flto -lm local/lib/libz.a local/lib/libdivsufsort64.a
#-fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wextra -Wall -fPIC -std=gnu99 -O2 -Ilocal/include/ -Isrc/ -g

EXEC = hera-T

SRC = src/alignment.c 				\
      src/barcode.c 				\
      src/bwt.c 				\
      src/dqueue.c 				\
      src/dynamic_alignment.c 			\
      src/genome.c 				\
      src/get_buffer.c 				\
      src/hash_table.c 				\
      src/index.c 				\
      src/io_utils.c 				\
      src/kmhash.c 				\
      src/opt.c 				\
      src/pthread_barrier.c     		\
      src/semaphore_wrapper.c 			\
      src/single_cell.c 			\
      src/utils.c 				\
      src/verbose.c 				\
      src/library_type.c 			\
      src/main.c

all:
	$(CC) -o $(EXEC) $(CFLAGS) $(SRC) $(LIBS)

clean:
	rm -f $(EXEC)
