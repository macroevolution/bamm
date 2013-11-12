all: bamm
	mkdir build
	mv src/bamm build/bamm

bamm:
	cd src && $(MAKE) $@

clean:
	cd src && $(MAKE) clean
	rm -rf build
