![ripta_logo](./ripta_logo.png)

# Re-evaluating Introgression site Patterns Through Ancestral alleles

**NOTE:** This is a beta version of `RIPTA` and is subject to change, all equations implemented are described in `appendix/RIPTA-equations.pdf` and the weighted-block Jackknife procedure is described in `appendix/RIPTA-weighted-block-jackknife-procedure.pdf`. If you have any questions or feature requests, please leave a detailed issue. If you use `RIPTA` in your work, please cite doi: https://doi.org/10.1101/2022.12.02.518851.

---

## Installation

These steps only need to be performed once.

### 1. Create the virtual environment

```bash
python -m venv RIPTA
source RIPTA/bin/activate
```

### 2. Install all dependencies

```bash
python -m pip install -r requirements.txt
```

---

## Usage

### Command line

```bash
source RIPTA/bin/activate
ripta --yaml_path config.yaml
```

### Python

```python
from ripta import read_config, process_vcfs, compile_results

config = read_config("config.yaml")
blocks = process_vcfs(config)
compile_results(config, blocks)
```