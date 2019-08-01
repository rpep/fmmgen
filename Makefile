install:
	pip3 install -e . --user

runtests:
	cd tests && py.test -v

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
