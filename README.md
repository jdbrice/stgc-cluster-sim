# stgc-cluster-sim 
---

## Get the code
```bash
git clone https://github.com/jdbrice/stgc-cluster-sim.git
cd stgc-cluster-sim
```

## build - must have ROOT on path
```bash
make 
```

## Run
```bash
./sim
```

## Analyze
Produces `output.root` that contains:
- `hGEN0` - true positions of the center of clusters
- `hGEN1` - clusters with realistic sigma
- `hGEN2` - clusters with realistic sigma + noise
- `hDIG`  - clusters with realistic sigma + noise digitized (to XY strip level) 