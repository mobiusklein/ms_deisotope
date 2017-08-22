dev:
	python setup.py develop

test:
	py.test -v  ms_deisotope --cov=ms_deisotope --cov-report=html --cov-report term

retest:
	py.test -v ms_deisotope --lf