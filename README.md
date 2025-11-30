# Comms_LAB — SDR Experiments & Documentation

This repository contains a set of software-defined radio (SDR) experiments, code, and documentation developed for the Communications Laboratory (IITH). The work focuses on OFDM-based communications using USRP and Adalm Pluto platforms and investigates physical-layer security (Alice–Bob–Eve scenarios).

**Compiled manual:** `SDR_MANUAL.pdf` is uploaded in this folder — open that PDF to view the full manual and hands-on guides.

This README gives a bird's-eye view of the project, where to find code and data, and how to build the project documentation (PDF).

---

## Quick summary
- Purpose: Implement and evaluate OFDM and MIMO PHY systems on heterogeneous SDR platforms, and study physical-layer security countermeasures.
- Primary devices used: Ettus USRP (N200/N210 family) and Analog Devices Adalm Pluto.
- Languages / tools used: MATLAB (primary), GNU Radio / Python (optional), LaTeX for documentation.

---

## Top-level layout (where to find what)

<!-- 
- `main.tex` — Main session document / template (session documentation, TOC, and layout). -->
- `cross_platform_pls.tex` — Detailed, standalone LaTeX report describing cross-platform physical-layer security and USRP↔Pluto interoperability (created by the tooling).
- `TRY_2_USRP/`, `Try2_mod/`, `Try2_TWOTX/`, `USRP_pluto/` — Project folders with MATLAB code, saved data, and hardware interface scripts. Each folder typically contains:
  - `Global_Parameters_PLS.m` — centralized experiment configuration (sampling rate, FFT size, gains, preamble definitions).
  - `OFDM_TX_*.m`, `OFDM_RX_*.m` — transmitter and receiver PHY implementations (OFDM chain).
  - `Hardware_TX_*.m`, `Hardware_RX_*.m` — hardware-specific wrappers / initialization scripts for USRP or Pluto.
  - `TX_dual_pluto.m` / `TX_signal_*.mat` — example transmission scripts and saved IQ payloads.
  - `corr_code.m`, `oversamp.m`, `setstate0_TX.m`, `setstate0_RX.m` — helper functions for synchronization, shaping and device state management.
- `USRP_pluto/` — contains `TX_dual_pluto.m` and receiver TX/RX scripts specific to Pluto-USRP experiments and example payload data files.

---

## Overall documentation (PDF)
We maintain a compiled documentation PDF that covers the projects and provides detailed methodology and experimental results.
- `cross_platform_pls.pdf` — Cross-platform PLS report (generated from `cross_platform_pls.tex`).
- `main.pdf` or `session_documentation.pdf` — Compiled `main.tex` if you prefer the session-style layout.

If the PDF(s) are not present, build them locally using the LaTeX instructions below.

---

## How to build the documentation (Windows PowerShell)
Prerequisites:
- LaTeX distribution installed (MiKTeX or TeX Live recommended) with `pdflatex` available in PATH.
- If `tikz`/`tcolorbox` packages are missing, the distribution should prompt automatic install or you may install them manually.

Open PowerShell in the `DOCS` folder and run:

```powershell
cd 'path'
# Build the cross-platform PLS report
pdflatex cross_platform_pls.tex
pdflatex cross_platform_pls.tex

# Optionally build the main session document
pdflatex main.tex
pdflatex main.tex
```

Notes:
- Run `pdflatex` twice to ensure table of contents and cross-references are resolved.
- If an `output.pdf` or `cross_platform_pls.pdf` exists and you want a fresh build, rename or delete the old PDF before building.

---

## How to run experiments (high-level)
Prerequisites:
- MATLAB (with Instrument Control / Support Packages if you use MATLAB's USRP or Pluto interfaces) OR GNU Radio / UHD (for USRP) and libiio/ADI support for Pluto.
- USRP UHD drivers for Ettus devices (install UHD and test with `uhd_find_devices` / `uhd_usrp_probe`).
- libiio and Pluto support (or the MATLAB PlutoSDR support package) for Adalm Pluto.
- Appropriate RF licenses and safety measures when transmitting in over-the-air tests.

Typical workflow:
1. Edit `Global_Parameters_PLS.m` to set `fs`, `fc`, `Nfft`, gains and device flags.
2. Start receiver scripts first (e.g., `OFDM_RX_Bob.m`) to capture preamble and avoid missing the first transmitted frame.
3. Run the transmitter script (e.g., `OFDM_TX_Alice.m` or `TX_dual_pluto.m`).
4. Collect logs and post-process with included analysis scripts (BER, CIR, secrecy metrics).

Troubleshooting tips are included in `cross_platform_pls.tex` and the `main.tex` document; consult them for synchronization and buffering issues.

---

## Files of interest (quick map)
- `Global_Parameters_PLS.m` — central configuration (must be consistent across TX/RX)
- `OFDM_TX_alice.m`, `OFDM_RX_bob.m` (naming varies by folder) — PHY examples
- `Hardware_TX_*.m`, `Hardware_RX_*.m` — hardware abstraction for USRP/Pluto
- `TX_dual_pluto.m` — Pluto dual-transmit example (in `USRP_pluto/`)
- `*.mat` — recorded payloads, channel logs and preamble sequences

---

## Physical-layer security focus
This repo includes experiments designed to evaluate secrecy capacity, artificial-noise injection, beamforming strategies, and the practical gap between theory and hardware reality. The key research goal is to examine Alice–Bob–Eve scenarios under heterogeneous hardware conditions (e.g., USRP↔Pluto) and report reproducible results.



## Contact
Communications Laboratory — IITH
Email: ee23btech11015@iith.ac.in or ee23btech11032@iith.ac.in

Note: This project was originally based on an existing open-source implementation. We have modified, extended, and restructured the code for my learning and experimentation. Original source: https://github.com/MeowLucian/SDR_Matlab_OFDM_802.11n.
---


