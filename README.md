# HAPS Guidance and Survey Optimization

MATLAB framework for modeling and optimizing autonomous high-altitude platform station (HAPS) surveys of Antarctic grounding zones using synthetic aperture radar (SAR).

## Overview

This repository contains modeling and optimization tools developed to investigate the use of a solar-powered high-altitude platform station (HAPS) for repeated SAR observations of Antarctic grounding zones. The project focuses on designing autonomous flight trajectories that balance scientific observation requirements with the physical and operational constraints of a long-endurance aircraft.

The framework combines models of Antarctic ice-shelf flexure, SAR and InSAR measurement requirements, atmospheric conditions, aircraft energy consumption, trajectory feasibility, and survey scheduling. These components are ultimately integrated into an optimization framework that determines **where and when the aircraft should collect observations**.

The primary components of the project include:

* **Ice-shelf flexure modeling** to determine the spatial extent of tidally driven grounding-zone deformation.
* **Temporal sampling analysis** to establish revisit requirements for resolving tidal signals.
* **InSAR uncertainty modeling** to evaluate the expected quality of displacement measurements.
* **Atmospheric modeling** using ERA5 wind fields to characterize environmental conditions encountered by the aircraft.
* **Aircraft energy modeling** to simulate solar power generation, propulsion requirements, and battery state of charge.
* **Reachability analysis** to determine SAR swath range and aircraft velocity.
* **Trajectory generation** for constructing candidate SAR survey paths.
* **Survey optimization** to balance coverage, revisit time, measurement uncertainty, and aircraft constraints.

## Repository Structure

```text
HAPS_Guidance/
│
├── ERA5_visualization/
├── InSAR_uncertainty/
├── beam_bending/
├── energy/
├── functions/
├── legacy/
├── optimization/
├── reachability_analysis/
├── temporal_requirement/
└── trajectories/
```

### `ERA5_visualization/`

Tools for processing and visualizing ERA5 atmospheric data, with an emphasis on the high-altitude wind fields relevant to aircraft trajectory planning.

### `InSAR_uncertainty/`

Models and analyses used to characterize uncertainty in airborne InSAR measurements and investigate how observation geometry and aircraft state affect measurement quality.

### `beam_bending/`

Models of tidally induced ice-shelf flexure near Antarctic grounding lines. These calculations are used to estimate the spatial extent of measurable tidal deformation and help define the survey region of interest.

### `energy/`

Aircraft energy-balance models, including aerodynamic power requirements, solar power generation, and battery state-of-charge calculations.

### `functions/`

Shared MATLAB functions used by multiple components of the project.

### `legacy/`

Older scripts and previous implementations retained for reference but not intended to represent the current workflow.

### `optimization/`

Development of the survey optimization framework used to determine the timing and location of SAR observations while balancing scientific objectives and aircraft constraints.

### `reachability_analysis/`

Analysis of SAR swath dimensions and groundspeed

### `temporal_requirement/`

Analysis of the temporal sampling and revisit requirements needed to resolve tidally driven signals.

### `trajectories/`

Generation and analysis of candidate aircraft trajectories and SAR survey patterns.

## Dependencies and Data

This project is written in **MATLAB** and relies on several external packages and geophysical datasets. These are not included directly in the repository and must be installed or downloaded separately.

### MATLAB Packages

* [Antarctic Mapping Tools (AMT)](https://github.com/chadagreene/Antarctic-Mapping-Tools) — Antarctic mapping, coordinate transformations, and access to Antarctic geospatial datasets.
* [Solar Position Calculator](https://www.mathworks.com/matlabcentral/fileexchange/58405-solar-position-calculator) — calculation of solar zenith and azimuth angles for the aircraft energy model.

Some analyses also require the MATLAB **Statistics and Machine Learning Toolbox** and **Image Processing Toolbox**.

### Datasets

* [ERA5 Hourly Pressure-Level Reanalysis](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels) — atmospheric wind fields used for high-altitude wind and reachability analyses.
* [MEaSUREs BedMachine Antarctica](https://nsidc.org/data/nsidc-0756/versions/4) — ice thickness, bed elevation, surface elevation, and ice/ocean masks.
* [MEaSUREs InSAR-Based Antarctica Ice Velocity Map](https://nsidc.org/data/nsidc-0484/versions/2) — Antarctic surface ice velocity used in geographic and survey-region analyses.
* [CATS2008 Antarctic Tide Model](https://www.usap-dc.org/view/dataset/601235)— tidal amplitudes and phases used to characterize tidal forcing and temporal observation requirements.

Large datasets are stored locally and are not tracked in this repository. Local file paths may need to be updated before running individual scripts.

## Requirements

The project is written in **MATLAB**.

Several components rely on external MATLAB packages and geophysical datasets. These include tools and datasets for Antarctic mapping, ice velocity, ice-sheet geometry, atmospheric conditions, and tidal modeling.

Specific dependencies and data sources will be documented as the project develops.

## Research Context

This repository is being developed as part of research into autonomous high-altitude SAR observations of Antarctic grounding zones.

The broader objective is to investigate how persistent, solar-powered aircraft could complement satellite observations by enabling adaptive and high-frequency measurements of rapidly evolving regions of the Antarctic Ice Sheet.

## Author

**Jeremy Wang**
Dartmouth College

Research conducted in collaboration with the Minchew Research Group at the California Institute of Technology.
