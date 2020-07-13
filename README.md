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
- `hDIG{N}` - to XY strip level
    - N=0 : noise only
    - N=1 : cluster only
    - N=2 : noise + cluster
- `hDIG{N}L{iMod}` - for module iMod = 1, 2, 3, 4 (top right, top left, bottom left, bottom right)
- `hDIG{N}L{iMod}hG{iG}` - projected 1D measurements for horizontal strips level
    - N = 3 or 4 (3 is projection, 4 applies saturation)
    - iMod = module 1, 2, 3, 4 (see above)
    - iG = strip group 0, 1, 2
- `hDIG{N}L{iMod}vG{iG}` - projected 1D measurements for vertical strips level






### Notes for me:

Plots generated with:

```bash
plot1d.sh --f=output.root --n=hstgc --d=colz --stat=0 --logz=1 --ext=png --Plot.Histo:title=";x (mm); y (mm)"

plot1d.sh --f=output.root --n=hGEN1 --d=colz --stat=0 --logz=1 --ext=png --Plot.Histo:title=";x (mm); y (mm)"

plot1d.sh --f=output.root --n=hGEN2 --d=colz --stat=0 --logz=1 --ext=png --Plot.Histo:title=";x (mm); y (mm)"

plot1d.sh --f=output.root --n=hDIG0 --d=colz --stat=0 --logz=1 --ext=png

plot1d.sh --f=output.root --n=hDIG1 --d=colz --stat=0 --logz=1 --ext=png
```