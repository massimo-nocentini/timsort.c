
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

# https://depth-first.com/articles/2019/10/16/compiling-c-to-webassembly-and-running-it-without-emscripten/
wasm-32:
	clang --target=wasm32-unknown-wasi -Wl,--import-memory -nostartfiles -Wl,--import-memory -Wl,--export-all -Wl,--no-entry -L/usr/local/lib/wasm32-wasi -DNDEBUG -o web/timsort.wasm timsort.c 

# if you have spiffy installed, otherwise `sudo chicken-install spiffy`
web:
	csi -e "(import spiffy) (start-server)" 