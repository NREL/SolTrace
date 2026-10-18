![Logo](https://github.com/NLR-SolTrace/SolTrace/tree/develop/api/assets/logo.png)

---

# PySolTrace

Python bindings and utilities for the NLR SolTrace&trade; ray-tracing engine.

---

## Installation

### Install from PyPI

```bash
pip install pysoltrace
```

### Verify Installation

```python
import pysoltrace

print(pysoltrace.__version__)
```

---

## Requirements

- Python 3.10+

---

## Examples

See [SolTrace Examples repository](https://github.com/NLR-SolTrace/SolTrace-Examples) for examples using pysoltrace.

---

## Development

See the main [README](https://github.com/NLR-SolTrace/SolTrace/blob/develop/README.md) for building SolTrace. For pysoltrace, set the `SOLTRACE_BUILD_API` CMake option to `ON`.

Navigate to `SolTrace/api`. Create virtual environment and install dependencies:

```bash
python -m venv .venv
.venv/Scripts/activate
pip install -r requirements.txt
```

To run tests, with the virtual environment activated, navigate to `SolTrace/api` and run:

```bash
python ./testing.py
```

---

## Building from Source

Navigate to `SolTrace/api/pysoltrace` and run:

```bash
python -m build --wheel
```

---

## SolTrace-pysoltrace Version Compatibility

| pysoltrace | SolTrace | Python |
|------------|----------|--------|
| 0.1.0 | 4.0.0-beta_v2 | 3.10+ |

---

## Contributing

See the main [CONTRIBUTING.md](https://github.com/NLR-SolTrace/SolTrace/blob/develop/CONTRIBUTING.md) for how to contribute to pysoltrace.

---

## Citation

We appreciate your use of SolTrace, and ask that you appropriately cite the software in exchange for its open-source publication. Please use one of the following references in documentation that you provide on your work. For general usage citations, the preferred option is:

```text
Wendelin, T. (2003). "SolTRACE: A New Optical Modeling Tool for Concentrating Solar Optics." Proceedings of the ISEC 2003: International Solar Energy Conference, 15-18 March 2003, Kohala Coast, Hawaii. New York: American Society of Mechanical Engineers, pp. 253-260; NREL Report No. CP-550-32866.
```

For citations in work that involves substantial development or extension of the existing code, the preferred option is:
```text
Wendelin, T., Wagner, M.J. (2018). "SolTrace Open-Source Software Project: github.com/NREL/SolTrace". National Laboratory of the Rockies. Golden, Colorado.
```

---

## License

See [LICENSE.md](https://github.com/NLR-SolTrace/SolTrace/tree/develop/coretrace/stapi_v2/pysoltrace/LICENSE.m) for details.
