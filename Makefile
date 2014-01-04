bamm:
	mkdir -p build && cd build && cmake ../ && make -j8

clean:
	rm -rf build
