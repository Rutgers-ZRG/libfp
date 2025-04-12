# Contributing to libfp

Thank you for your interest in contributing to libfp! This document provides guidelines and instructions for contributing.

## Code of Conduct

Please note that this project adheres to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## How to Contribute

### Reporting Bugs

If you find a bug in the software, please create an issue on our [GitHub Issues](https://github.com/Rutgers-ZRG/libfp/issues) page with the following information:

- A clear, descriptive title
- A detailed description of the issue
- Steps to reproduce the bug
- Expected behavior and actual behavior
- Your environment (OS, Python version, etc.)

### Suggesting Enhancements

We welcome suggestions for enhancements to libfp. Please create an issue on our [GitHub Issues](https://github.com/Rutgers-ZRG/libfp/issues) page with the following information:

- A clear, descriptive title
- A detailed description of the proposed enhancement
- Any relevant examples or use cases

### Pull Requests

We welcome pull requests for bug fixes, enhancements, and documentation improvements. Here's how to create a pull request:

1. Fork the repository
2. Create a new branch for your changes
3. Make your changes in your fork
4. Run the tests to ensure your changes don't break existing functionality
5. Create a pull request to our main branch

Please include the following in your pull request:

- A clear, descriptive title
- A detailed description of the changes
- Any relevant issue numbers
- Tests for your changes

## Development Setup

To set up the development environment:

1. Clone the repository
    ```bash
    git clone https://github.com/Rutgers-ZRG/libfp.git
    cd libfp
    ```

2. Install the development dependencies
    ```bash
    pip install -e ".[dev]"
    ```

3. Run the tests
    ```bash
    pytest
    ```

## Testing

All new code should include appropriate tests. We use pytest for testing.

To run the tests, use:

```bash
pytest
```

## Documentation

Please document all public classes, methods, and functions using [NumPy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).

## License

By contributing to libfp, you agree that your contributions will be licensed under the project's [MIT License](LICENSE).