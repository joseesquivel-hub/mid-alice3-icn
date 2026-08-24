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
