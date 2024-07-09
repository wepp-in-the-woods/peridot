# peridot
Programmable Environmental Rust Interface for Drainage &amp; Operational Topography.

This library provides topological processing routines for __abstracting__ hillslopes and channels for WEPP from TOPAZ ([TOpographic PArameteriZation](https://github.com/rogerlew/topaz)).

## Subcatchment Delineation

WEPPcloud uses the subcatchments identified from TOPAZ (SUBWTA.ARC) to define rectangular hillslopes (length and width) with slope profiles. This is accomplished by walking the flowpaths in the subcatchment to identify a collection of flowpath slope profiles. 
Then the flowpaths are aggregated to define a singular slope profile for the subcatchment. For left and right subcatchments the width is defined as the length of the cooresponding channel that the subcatchments drain into. The length of the top subcatchments is defined from the flowpath lengths and the width is determined from the area and length.

## Channel Delineation

The channel slope profiles are defined by walking the channel and the width is determined from the channel's upland area ([#268](https://github.com/rogerlew/wepppy/issues/268)).

## Catchment Trace

This library also contains a catchment tracing algorithm that delineates a single hillslope from a pourpoint location (e.g. culvert). The algorithm identifies the uparea pixels and then walks down to the defined pourpoint. Once the uparea pixels are identified the subcatchment routine is used to define abstract rectangular hillslope profiles for WEPP.


## build instructions

```bash
cargo build --release
```


## wepppy integration

```bash
cp /workdir/peridot/target/release/abstract_watershed /workdir/wepppy/wepppy/topo/peridot/bin/
```
