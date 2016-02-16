release : build/release/msmc2

all: debug release decode

debug : build/debug/msmc2

all : debug release

unittest: build/debug/unittest
	build/debug/unittest

build/debug/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -debug -L-lgsl -L-lgslcblas -odbuild/test -of$@ $^

build/release/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -O -L-lgsl -L-lgslcblas -odbuild/test -of$@ $^

build/debug/unittest : model/*.d
	dmd -unittest -L-lgsl -L-lgslcblas -main -odbuild/debug -ofbuild/debug/unittest $^

decode : build/decode

build/decode : model/*.d decode.d 
	dmd -O -L-lgsl -L-lgslcblas -odbuild -of$@ $^

clean :
	rm -rf build/debug/* build/release/*

.PHONY : clean unittest

