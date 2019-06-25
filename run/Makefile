CC = gcc
CFLAGS = -O3 -fPIC -Wall
LDFLAGS = -shared

all: _get_data.so

_get_data.so: get_data.c
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $<

clean:
	rm -f _get_data.so
