dev:
	pip install -e . --no-build-isolation -v

clean:
	@-rm ms_deisotope/_c/*.pyd
	@-rm ms_deisotope/_c/*/*.pyd
	@-rm ms_deisotope/_c/*.pdb
	@-rm ms_deisotope/_c/*/*.pdb
	@-rm -rf build

test:
	py.test -v  ms_deisotope --cov=ms_deisotope --cov-report=html --cov-report term

retest:
	py.test -v ms_deisotope --lf


update-cv-lists:
	python -m cogapp -r ms_deisotope/data_source/metadata/software.py \
		ms_deisotope/data_source/metadata/activation.py \
		ms_deisotope/data_source/metadata/data_transformation.py \
		ms_deisotope/data_source/metadata/instrument_components.py \
		ms_deisotope/data_source/metadata/file_information.py \
		ms_deisotope/data_source/metadata/scan_traits.py
	python -m autopep8 -i --max-line-length 80 ms_deisotope/data_source/metadata/software.py \
		ms_deisotope/data_source/metadata/activation.py \
		ms_deisotope/data_source/metadata/data_transformation.py \
		ms_deisotope/data_source/metadata/instrument_components.py \
		ms_deisotope/data_source/metadata/file_information.py \
		ms_deisotope/data_source/metadata/scan_traits.py


update-docs:
	git checkout gh-pages
	git pull origin master
	cd docs && make clean html
	git add docs/_build/html -f
	git commit -m "update docs"
	git push origin gh-pages
	git checkout master