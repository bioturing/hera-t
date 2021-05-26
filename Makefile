CC = gcc

CPP = cpp

AR = ar

LIBS = -pthread -flto -lm local_unix/lib/libz.a local_unix/lib/libdivsufsort64.a
#-fsanitize=undefined,address

CFLAGS = -Wfatal-errors -Wextra -Wall -fPIC -std=gnu99 -O2 -Ilocal_unix/include/ -Isrc/
#-g

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

OBJ = $(SRC:.c=.o)

DEP = $(OBJ:.o=.d)

.PHONY: all
all: $(EXEC)

.PHONY: debug
debug: LIBS += -fsanitize=undefined,address
debug: CFLAGS += -g
debug: cleanall
debug: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

-include $(DEP)

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -rf $(OBJ) $(EXEC)

.PHONY: cleandep
cleandep:
	rm -f $(DEP)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ) $(EXEC) $(DEP)

