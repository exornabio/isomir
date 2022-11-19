CC := g++
CFLAGS := -Wall -g  -I../include
TARGET := isomir
object := main.o config.o detection.o levenshtein.o

.PHONY: clean

$(TARGET): $(object)
	$(CC) $^ -o $@ -static -g -lboost_system -lboost_filesystem
	cp $@ ../script

%.o: %.cpp
	$(CC) $< -o $@ -c -std=c++11 -O3 

clean :
	rm *.o isomir

