
linux:
	gcc -O2 -std=gnu99 -fPIC -Wall -g -I/usr/local/include/ -s -shared -o libtimsort.so timsort.c

install:
	cp timsort.h /usr/local/include/
	mv libtimsort.so /usr/local/lib/
