test:
	py.test -v  ms_deisotope --cov=ms_deisotope --cov-report=html

retest:
	py.test -v ms_deisotope --lf