CXX_COMPILER = g++
CXX_FLAGS =

src_dir = ~/projects/bamm/src/
gtest_include_dir = ~/gtest-1.7.0/include/
gtest_src_dir = ~/gtest-1.7.0/src/

test_files = *.cpp
src_files = $(src_dir)/[A-Z]*.cpp    # Excludes main.cpp

test-all: $(test_files) $(src_files)
	$(CXX_COMPILER) $(CXX_FLAGS) -I $(src_dir) -I $(gtest_include_dir) \
		-o test-all $(test_files) $(src_files) \
		$(gtest_src_dir)/gtest-all.o $(gtest_src_dir)/gtest_main.o

clean:
	rm test-all
