# SCS module containers

This module uses a single pre-built container image:

- `ghcr.io/katwre/scs:latest`

The Dockerfile clones upstream SCS (`https://github.com/chenhcs/SCS.git`) and creates the conda environment from `/app/SCS/environment.yml`.
To keep builds reproducible, the Dockerfile pins SCS to commit `44bc89b2fccdc51cbd59ecfd50d1ff5a2134f378` by default.

## Build examples

Build image:

```bash
docker build -f modules/local/scs/Dockerfile -t ghcr.io/katwre/scs:latest .
```

Override pinned commit if needed:

```bash
docker build -f modules/local/scs/Dockerfile --build-arg SCS_COMMIT=<commit-sha> -t ghcr.io/katwre/scs:latest .
```

Push:

```bash
docker push ghcr.io/katwre/scs:latest
```
