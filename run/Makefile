CC = gcc
CFLAGS = -O3 -fPIC -Wall
CXX = g++
CXXFLAGS = -O3 -fPIC -Wall
LDFLAGS = -shared

all: _get_data.so _summary_plot_averages.so

_get_data.so: get_data.c
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $<

_summary_plot_averages.so: summary_plot_averages.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $<

clean:
	rm -f _get_data.so _summary_plot_averages.so
