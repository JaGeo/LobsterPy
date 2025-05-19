# Developer Installation

Install LobsterPy from source, by cloning the repository via [github](https://github.com/JaGeo/LobsterPy.git)

```bash
git clone https://github.com/JaGeo/LobsterPy.git
cd LobsterPy
pip install -e .[featurizer,docs,tests,dev]
```
This will install LobsterPy will all dependencies for tests, pre-commit and docs building.

## Running unit tests

Unit tests can be run from the source folder using `pytest`. 

```bash
pytest
```
This will run all the tests.

To get a detailed report of test coverage you can use following command
```bash
pytest --cov=lobsterpy --cov-report term-missing --cov-append
```

If you feel test execution takes too long locally, you can speedup the execution using [pytest-xdist](https://pypi.org/project/pytest-xdist/). Install this in library in your environment using

```bash
pip install pytest-xdist
```

Once installed, you can now use multiple processors to run your tests. For example, if you want to use 8 processors to run tests in parallel, run

```bash
pytest -n 8
```

We rely on pytest-split to run tests in parallel on github workflow, thus it is necessary to update the test-durations files in the repository, incase you add new tests. To generate this file, use

```bash
pytest --cov=lobsterpy --cov-append --splits 1 --group 1 --durations-path ./tests/test_data/.pytest-split-durations --store-durations
```

## Building the documentation locally

The atomate2 documentation can be built using the sphinx package.

The docs can be built to the `_build` directory:

```bash
sphinx-build -W docs _build
```

