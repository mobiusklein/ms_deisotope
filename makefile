dev:
	pip install -e . --no-build-isolation -v

clean:
	@-rm src/ms_deisotope/_c/*.pyd
	@-rm src/ms_deisotope/_c/*/*.pyd
	@-rm src/ms_deisotope/_c/*.pdb
	@-rm src/ms_deisotope/_c/*/*.pdb
	@-rm -rf build

test:
	py.test -v ./tests --cov=ms_deisotope --cov-report=html --cov-report term

retest:
	py.test -v ./tests --lf --pdb


update-cv-lists:
	python -m cogapp -r src/ms_deisotope/data_source/metadata/software.py \
		src/ms_deisotope/data_source/metadata/activation.py \
		src/ms_deisotope/data_source/metadata/data_transformation.py \
		src/ms_deisotope/data_source/metadata/instrument_components.py \
		src/ms_deisotope/data_source/metadata/file_information.py \
		src/ms_deisotope/data_source/metadata/scan_traits.py
	python -m autopep8 -i --max-line-length 80 src/ms_deisotope/data_source/metadata/software.py \
		src/ms_deisotope/data_source/metadata/activation.py \
		src/ms_deisotope/data_source/metadata/data_transformation.py \
		src/ms_deisotope/data_source/metadata/instrument_components.py \
		src/ms_deisotope/data_source/metadata/file_information.py \
		src/ms_deisotope/data_source/metadata/scan_traits.py


update-docs:
	git checkout gh-pages
	git pull origin master
	cd docs && make clean html
	git add docs/_build/html -f
	git commit -m "update docs"
	git push origin gh-pages
	git checkout master