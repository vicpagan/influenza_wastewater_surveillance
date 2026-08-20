# the compiler: gcc for C program
CC = gcc

# compiler flags:
# -Wall turns on most, but not all, compiler warnings
# -O3 is for compiler optimizations
# -g adds debugging information to the executable file
CFLAGS = -Wall -O3
DBGCFLAGS = -Wall -O0 -g -fsanitize=address

# -lz links the gz library, -lpthread links pthreads
LIBS = -lz -lpthread -std=gnu99

# sources
SOURCES = eliminate_strains.c calculate_proportions.c calculate_allele_freq.c build_mismatch_matrix.c options.c msa.c sam.c bowtie_alignment.c clean_reads.c file_utils.c align_reference.c
EXTERNAL = external/needleman_wunsch.c external/alignment.c external/alignment_scoring.c external/covid_em.c external/hashmap.c
OBJECTS = $(SOURCES:.c=.o)

# the build target executable:
TARGET = eliminate_strains

all: $(TARGET)

$(TARGET): $(SOURCES) $(EXTERNAL)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES) $(EXTERNAL) $(LIBS)

debug: $(SOURCES) $(EXTERNAL)
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(SOURCES) $(EXTERNAL) $(LIBS)

clean:
	$(RM) $(TARGET)