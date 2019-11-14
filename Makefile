all:
	make install
	make test
	make ipython

docker-all:
  # docker-build also installs software in container
	make docker-build
	make docker-test
	@echo "Build, install and tests run. Looking good."


install:
	python3 -m pip install -e .

test:
	cd tests && py.test -v

ipython:
	ipython

# To use Docker container for building and testing
# build docker image locally, needs to be done first

docker-build:
	docker build -t fmmgen -f Dockerfile .

docker-test:
	docker run --rm -v `pwd`:/io fmmgen make test

# for interactive exploration
docker-ipython:
	docker run --rm -ti -v `pwd`:/io fmmgen ipython

# for diagnostic purposes and convenience
docker-bash:
	docker run --rm -ti -v `pwd`:/io fmmgen bash


clean:
	rm -rf tests/operators_decl.pxd
	rm -rf tests/operators_wrap.pyx
	rm -rf tests/operators_wrap.pyxbld
	rm -rf tests/operators.h
	rm -rf tests/operators.c
	rm -rf tests/__pycache__
	rm -rf *.egg-info
	rm -rf __pycache__
	rm -rf *~
	rm -rf tests/LinearDipole.c
	rm -rf tests/LinearDipole.h
	rm -rf tests/LinearQuadrupole.c
	rm -rf tests/LinearQuadrupole.h
	rm -rf tests/MonopoleOrigin.c
	rm -rf tests/MonopoleOrigin.h
	rm -rf tests/QuadrupoleTwoDipoles.c
	rm -rf tests/QuadrupoleTwoDipoles.h
