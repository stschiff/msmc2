all: debug release

debug : build/debug/msmc2

release : build/release/msmc2

unittest: build/debug/unittest
	build/debug/unittest

build/debug/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -debug -odbuild/test -of$@ $^

build/release/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -O -release -odbuild/test -of$@ $^

build/debug/unittest : model/*.d
	dmd -unittest -main -odbuild/debug -ofbuild/debug/unittest $^
	
clean :
	rm -rf build/debug/* build/release/*

.PHONY : clean unittest

