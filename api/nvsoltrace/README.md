![Logo](https://github.com/NLR-SolTrace/SolTrace/tree/develop/api/assets/logo.png)

---

# NVSolTrace

NVIDIA Optix plugin for PySolTrace, the python bindings for the NLR SolTrace&trade; ray-tracing engine. Must have an NVIDIA GPU to run this code. This just packages the GPU related code seperately from the CPU code, and it will be automatically included if installed.

---

## Installation

### Install from PyPI

```bash
pip install pysoltrace[optix]
```

### Verify Installation

```python
import nvsoltrace

print(nvsoltrace.__version__)
```

---

## Requirements

- Python 3.10+

---

## Examples

See [SolTrace Examples repository](https://github.com/NLR-SolTrace/SolTrace-Examples) for examples using the OPTIX runner from pysoltrace.

---

## Development

See the main [README](https://github.com/NLR-SolTrace/SolTrace/blob/develop/README.md) for building SolTrace. For nvsoltrace, set both `SOLTRACE_BUILD_API` and `SOLTRACE_BUILD_OPTIX_SUPPORT` CMake options to `ON`.

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

If nvsoltrace is successfully found, the OptixRunner tests will run.

---

## Building from Source

Install build tool. Navigate to `SolTrace/api/nvsoltrace` and run:

```bash
python -m pip install build
python -m build --wheel
```

---

## Version Compatibility

| nvsoltrace | pysoltrace | SolTrace | Python |
|------------|------------|----------|--------|
| 0.0.1 | 0.0.1 | 4.0.0-beta_v2 | 3.10+ |

---

## Contributing

See the main [CONTRIBUTING.md](https://github.com/NLR-SolTrace/SolTrace/blob/develop/CONTRIBUTING.md) for how to contribute to nvsoltrace.

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

See [LICENSE.md](https://github.com/NLR-SolTrace/SolTrace/tree/develop/api/nvsoltrace/LICENSE.md) for details.
