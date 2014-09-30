msmc2 : build/msmc2

build/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -O -odbuild -ofbuild/msmc2 $^

build/test/msmc2 : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc2.d logger.d
	dmd -O -odbuild/test -ofbuild/test/msmc2 $^

clean :
	find build -type f -delete

.PHONY : clean

