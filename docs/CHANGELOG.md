# CHANGELOG — MID Geometry ALICE 3 / AliceO2

## August 2026 — ICN-UNAM

---

### Commit 595090e — ALICE3Field.C correction
- Fixed B2 equation: B1 moved from denominator to numerator
  - Before: `B2 = -Rc*Rc / ((R2*R2 - R1*R1) * B1)` → -0.479 T
  - After:  `B2 = -B1*Rc*Rc / (R2*R2 - R1*R1)` → -1.619 T
- Rc set to 170 cm (R_out_coil per Ian DetectorConstruction.cc)
- Confirmed by Dr. Antonio Ortiz, August 2026
- Field activated via ALICE3_SIM_FIELD=ON + ALICE3_MAGFIELD_MACRO

---

### Commit 84b53c5 — Layer 1 half-length correction
- Layer 1 half-length corrected from 520 cm (arbitrary) to 525 cm
- Source: Ian DetectorConstruction.cc midLength_2 = 1.05 * absorberLength
- Confirmed by Scoping Document Table 16: total length 10.5 m
- Module step updated dynamically: 52.5 cm (525/10 modules)
- Verified with InspectShape: total length 525.00 cm

---

### Commit 9fe6ea6 — Absorber geometry (Antonio request)
- Absorber outer radius fixed to Rmax=290 cm constant in all sections
- Stepped Rmin: 220 cm central (|z|<300 cm), 245 cm external (|z|>300 cm)
- Verified with InspectShape: Nz=6, Rmax=290 constant in all planes
- Source: Ian DetectorConstruction.cc + Scoping Document Table 16

---

### Commit b3174b1 — Hit production fix
- Fixed sensor names for VMC resolution
- Simple sensor names allow correct hit registration per layer
- Verified: hits correctly distributed between ID0 and ID1

---

### Commit a03c7f0 — ProcessHits layer index fix
- Fixed physLay extraction from sensor name
- Layer index correctly extracted from MIDSensor_L{n}_... name pattern

---

### Commit 2d1f6f0 — Initial implementation
- IanMagnet geometry: Al walls, NbTi/Cu/Al coil, Al supports
- IanAbsorber: stepped geometry with constant outer radius
- ICNStepped MID layout: 2 layers, 16 staves each
- 0 overlaps verified with Nicolo's macro

---

## Key design decisions

| Parameter | Value | Source |
|-----------|-------|--------|
| Absorber Rmax | 290 cm | Ian code + Scoping Doc Table 16 |
| Absorber Rmin central | 220 cm | Ian code |
| Absorber Rmin external | 245 cm | Ian code (45 cm thickness) |
| Layer 0 radius | 301 cm | Scoping Doc Table 16 |
| Layer 1 radius | 311 cm | Scoping Doc Table 16 |
| Layer 0 half-length | 499 cm | Ian code |
| Layer 1 half-length | 525 cm | Ian code + Scoping Doc Table 16 |
| Staves per layer | 16 | Scoping Doc Table 16 |
| Bar width | 2.5 cm | Ian DetectorConstruction.cc |
| Bar spacing (φ) | 2 mm | Ian DetectorConstruction.cc |
| Bar spacing (z) | 1 mm | Ian DetectorConstruction.cc |
| Magnet coil Rin | 160 cm | Ian DetectorConstruction.cc |
| Magnet coil Rout | 170 cm | Ian DetectorConstruction.cc |
| Rc (field boundary) | 170 cm | R_out_coil, confirmed A. Ortiz |
| B1 | +2 T | Design spec |
| B2 | -1.619 T | Flux conservation, corrected Aug 2026 |

---

## Open questions (pending confirmation from A. Ortiz)

- Bar spacing: 2 mm implemented (from Ian), Antonio mentioned 1 mm — to confirm
- Mechanical support structure not yet implemented (bars floating)
- Filling factor correction (~2-8%) not yet modeled in simulation

---

## August 26, 2026 — Track Extrapolation Analysis

### HitMap_O2.root
- Generated from scratch using our validated geometry (A. Ortiz, Aug 2026)
- 6880 bars: L0=3680 (16x10x23), L1=3200 (16x10x20)
- Bar positions calculated analytically using Ian's polar formula with our radii
- Bar pitch = 5.2 cm both layers; L0 sensor length = 49.9 cm, L1 = 61.75 cm

### extrapolation.C
- Relativistic helix extrapolation from L0 (R=301 cm) to L1 (R=311 cm)
- Field B2 = -1.619 T in absorber; straight line between layers (no field)
- Results: sigma~1.8 cm (X,Y,Z), no systematic bias, 94.4% in +/-5 cm window
- Consistent with multiple scattering in 70 cm Fe absorber

### extrapolation_digitized_v4.C
- Digitization of L0 bar (X,Y center) and L1 bar (Z center)
- Search algorithm: phi filter -> Z range filter -> minimum XY distance (L0)
- Search algorithm: nearest stave in phi -> nearest bar in Z (L1)
- Results: 92.6% in +/-5 cm window — digitization costs 1.8 percentage points
- No bias after digitizing: means compatible with zero at <1.2 sigma

### plot_residuos_v4.C
- Final figure: blue=MC truth, red=digitized, green dashed lines at +/-5 cm
- RMS and mean reported for all 6 distributions

### Open questions / Next steps (Route B)
- Model real Z uncertainty of L0 starting point (sigma~29 cm from bar length)
- Current result uses MC truth Z as extrapolation starting point
- ΔZ RMS slightly below prediction (1.22 vs 1.50 cm): attributed to
  correlation between digitization error and residual — not a bug
- Run with pp-like generator (currently particle gun only)

---

## September 1, 2026 — Acceptance Analysis

### Fix: BoxGun.number=1
- Antonio identified that BoxGun.number was not set, defaulting to 10 muons/event
- Correct command requires BoxGun.number=1 for single muon per event
- Impact: hits/event from ~24 to 2.29 (consistent with prototype article: 2.7)

### Correct simulation command:
o2-sim -n 1000 --detectorList ALICE3 -g boxgen -m MI3 -m A3MAG -m A3ABSO \
  --configKeyValues "align-geom.mDetectors=none;Alice3PassiveBase.mLayout=3;\
Alice3PassiveBase.mDetLayout=2;MIDBase.mLayout=0;BoxGun.pdg=13;BoxGun.number=1;\
BoxGun.eta[0]=-1.2;BoxGun.eta[1]=1.2;BoxGun.prange[0]=1.5;BoxGun.prange[1]=5"

### Acceptance results (1000 events, g=1mm, 1.5-5 GeV, |eta|<1.2):
- L0 OR L1:              99.2%
- L0 AND L1:             97.3%
- L0 AND L1, same stave: 96.8% — consistent with Antonio's epsilon^2=0.96
