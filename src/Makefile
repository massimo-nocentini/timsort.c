
linux:
	gcc -O2 -g -std=c11 -fPIC -Wall -shared -o libtimsort.so timsort.c

mingw:
	gcc -O2 -g -std=c11 -fPIC -Wall -I/usr/local/include/ -shared -o libtimsort.dll timsort.c

macos:
	gcc -O2 -g -std=c11 -fPIC -Wall -dynamiclib -o libtimsort.dylib timsort.c

install:
	cp timsort.h /usr/local/include/
	mv libtimsort.so /usr/local/lib/

install-mingw:
	cp timsort.h /usr/local/include/
	mv libtimsort.dll /usr/local/lib/

install-macos:
	cp timsort.h /usr/local/include/
	mv libtimsort.dylib /usr/local/lib/
