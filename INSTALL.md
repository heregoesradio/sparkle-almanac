## Quickstart
### 1. Clone repository
```
git clone --recurse-submodules https://github.com/heregoesradio/sparkle-almanac.git
```

### 2. Install `heregoes` Conda environment
##### Intel (MKL)
```
conda env create -f heregoes/release/heregoes-env-intel.yml
```

##### AMD, ARM64 (OpenBLAS)
```
conda env create -f heregoes/release/heregoes-env-other.yml
```

### 3. Activate
```
conda activate heregoes-env
```

### 4. Import in Python
```python
from fast import FastSparkleAlmanac
from slow import SparkleAlmanac
```