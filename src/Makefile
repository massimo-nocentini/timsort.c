
linux:
	cc -O2 -std=c11 -fPIC -Wall -I/usr/local/include/ -s -shared -o libtimsort.so timsort.c

macos:
	cc -O2 -std=c11 -fPIC -Wall -dynamiclib -o libtimsort.dylib timsort.c

install:
	cp timsort.h /usr/local/include/
	mv libtimsort.so /usr/local/lib/

install-macos:
	cp timsort.h /usr/local/include/
	mv libtimsort.dylib /usr/local/lib/
