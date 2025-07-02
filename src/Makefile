
linux:
	clang -O3 -g -fPIC -Wall -shared -DNDEBUG -o libtimsort.so timsort.c

mingw:
	clang -O3 -g -fPIC -Wall -I/usr/local/include/ -shared -DNDEBUG -o libtimsort.dll timsort.c

macos:
	clang -O3 -g -fPIC -Wall -dynamiclib -DNDEBUG -o libtimsort.dylib timsort.c

install:
	cp timsort.h /usr/local/include/
	mv libtimsort.so /usr/local/lib/

install-mingw:
	cp timsort.h /usr/local/include/
	mv libtimsort.dll /usr/local/lib/

install-macos:
	cp timsort.h /usr/local/include/
	mv libtimsort.dylib /usr/local/lib/
