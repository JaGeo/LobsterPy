# Contributing to LobsterPy

We would love your input to further extend and improve our code. Contribution of any kind are welcome, for example it can be:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing or implementing new features
- Or simply updating the documentation

## Reporting bugs, getting help, and discussion

You can submit questions and bugs to the
[GitHub issues page](https://github.com/JaGeo/LobsterPy/issues).

When creating a bug report, its best if you include some of the pointers mentioned below:

- A quick summary and/or background.
- Steps to reproduce - be specific! **Provide sample code.**
- What was expected and what actually happened.
- Screenshots of any errors you encounter.

## Contributing code improvements or additions through Github

- Fork the repo and create your branch from master.
- If working using IDE, clone the repo and install the package in development mode `pip install -e lobsterpy[featurizer,dev,tests]` in a separate conda or virtual environment.
- run `pre-commit install`  (This will install all the hooks from `.pre-commit-config.yaml` and make adapting your code to match linting standards)
- Commit your improvements to your branch and push to your Github fork (repo).
- Make sure you write tests and update documentation when needed. We use pytest framework and check out our existing tests for inspiration.
- Try to reuse existing test data when possible.
- When you're finished, go to your fork and make a Pull Request (PR). It will
  automatically update if you need to make further changes.
- If you have raised a PR but are actively working on it, please add `[WIP] in the title.(This will let us know it is still not ready for review or merged)
