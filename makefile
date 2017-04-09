test:
	nosetests --with-coverage --with-timer --cover-package=ms_deisotope --cover-html --cover-html-dir=test_reports\
			  --logging-level=DEBUG -v --with-id ms_deisotope/test/

retest:
	nosetests --cover-package=ms_deisotope --logging-level=DEBUG -v --with-id --failed ms_deisotope/test/