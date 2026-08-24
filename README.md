# MID Geometry — ALICE 3 / AliceO2

Implementation of the Muon Identifier (MID) superconducting magnet geometry for the ALICE 3 upgrade, ported to the AliceO2 O2 framework.

**Status:** Verified — August 2026.
**Branch:** [mid-magnet-geometry](https://github.com/joseesquivel-hub/AliceO2/tree/mid-magnet-geometry)

---

## What this implements

- Superconducting magnet geometry: inner/outer Al walls, NbTi/Cu/Al coil, Al support structure
- Iron absorber: constant Rmax=290 cm, stepped Rmin (220 cm central, 245 cm external)
- Two MID detector layers: L0 at R=301 cm (|z|<499 cm), L1 at R=311 cm (|z|<525 cm), 16 flat staves each
- ALICE 3 magnetic field correction: B2 = -B1·Rc²/(R2²-R1²), Rc=170 cm, B2=-1.619 T

## Modified files (vs upstream/dev)

| File | Description |
|------|-------------|
| `src/Absorber.cxx` | Stepped absorber geometry |
| `src/Magnet.cxx` | Superconducting magnet layers |
| `src/PassiveBaseParam.h` | Absorber and magnet parameters |
| `src/Detector.cxx` | MID hit production, layer 1 half-length 525 cm |
| `src/MIDLayer.cxx` | Stave and module construction |
| `src/Detector.h` | MID detector header |
| `src/MIDLayer.h` | MIDLayer class with nBars, zOffset, staveWidth |
| `src/MI3BaseParam.h` | MID base parameters |
| `src/ALICE3Field.C` | ALICE 3 magnetic field macro (corrected) |

## Key results

- 0 overlaps verified with Nicolo's macro over full ALICE 3 geometry
- 24,048 hits / 1000 muon events (GEANT4, 5-50 GeV, |η|<1.2)
- Geometric acceptance: 80.4% (uniform generation in |η|<1.65)
- Flat plateau at 100% for |η|<1.15, sharp edge at |η|=1.28/1.30

## How to run

### Prerequisites
- AliceO2 environment loaded: `alienv enter O2/latest-o2`
- Copy files from `src/` to the corresponding paths in your AliceO2 tree
- Compile: `ninja install` from your O2 build directory

### Activate ALICE 3 field
```bash
export ALICE3_SIM_FIELD=ON
export ALICE3_MAGFIELD_MACRO=/path/to/src/ALICE3Field.C
```

### Run simulation
```bash
o2-sim -n 1000 --detectorList ALICE3 -g boxgen -m MI3 -m A3MAG -m A3ABSO \
  --configKeyValues "align-geom.mDetectors=none;Alice3PassiveBase.mLayout=3;\
Alice3PassiveBase.mDetLayout=2;MIDBase.mLayout=0;BoxGun.pdg=13;\
BoxGun.eta[0]=-1.2;BoxGun.eta[1]=1.2;BoxGun.prange[0]=5;BoxGun.prange[1]=50"
```

### Verify geometry
```bash
root -l -b -q -e 'TGeoManager::Import("o2sim_geometry.root"); \
  gGeoManager->GetVolume("Absorber")->InspectShape();'
```
Expected: Nz=6, Rmax=290 constant in all planes.

### Verify field
```bash
grep "ALICE3 magnetic field" o2sim_workerlog0
```
Expected: `Initializing ALICE3 magnetic field`

## Macros

| Macro | Description | Usage |
|-------|-------------|-------|
| `macros/pdg_count.C` | Particle type catalog in MID hits | `root -l -b -q pdg_count.C` |
| `macros/color_mid_geometry.C` | 3D colored geometry for JSROOT | `root -l -b -q color_mid_geometry.C` |

JSROOT viewer: https://root.cern/js/latest/

## Authors

Francisco Esquivel — ICN-UNAM / ALICE-CERN
Based on GEANT4 standalone simulation (ICN-UNAM)
