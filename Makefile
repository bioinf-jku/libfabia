CC=gcc -std=gnu99
CFLAGS=-O3 -march=native -Wall -fopenmp -DNDEBUG -I./ext/include/ -fPIC
CXXFLAGS=$(CFLAGS)
LDFLAGS=-O3 -march=native -flto -fopenmp
LIBS=-lhdf5 -lhdf5_hl -lm

# Tweak this if you have MKL available!
LIBS+=-L/usr/lib/lapack/ -llapack -lblas


# Tune these include paths to match your HDF5 installation
CFLAGS += -I/usr/include/hdf5/serial
LDFLAGS+= -L/usr/lib/x86_64-linux-gnu/hdf5/serial/


SOURCES=main.c fabia.c hdf5tools.c util.c fabia_approx.c
OBJECTS = $(patsubst %.c,%.o,$(SOURCES))
EXECUTABLE=fabia


all: $(SOURCES) $(EXECUTABLE) libfabia.so

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)

libfabia.so: fabia.o util.o fabia_approx.o interfaces.o updateui.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS) -shared

.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o myfabia.so fabia libfabia.so afabia
