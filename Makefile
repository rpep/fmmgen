all:
	make install
	make test
	make ipython

docker-all:
	make docker-build
	make docker-install
	make docker-test
	make docker-ipython
	make docker-bash


install:
	python3 -m pip install .

test:
	cd tests && py.test -v

ipython:
	ipython

# To use Docker container for building and testing
# build docker image locally, needs to be done first

docker-build:
	docker build -t fmmgen -f Dockerfile .

docker-install:
	docker run -v `pwd`:/io fmmgen make install

docker-test:
	docker run -v `pwd`:/io fmmgen make test

docker-ipython:
	docker run -ti -v `pwd`:/io fmmgen ipython

docker-bash:
	docker run -ti -v `pwd`:/io fmmgen bash


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
