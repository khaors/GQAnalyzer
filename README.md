# __GQAnalyzer__: R package to analyze Groundwater, Geothermal water and Gas Geochemistry information
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rves)](https://cran.r-project.org/package=GQAnalyzer)
[![HitCount](http://hits.dwyl.io/khaors/rves.svg)](http://hits.dwyl.io/khaors/GQAnalyzer)

The goal of GQAnalyzer is to be a tool to support the analysis of the chemical composition of different types of geological fluids (groundwater, geothermal waters and geothermal gases). This package includes several functions to create plots routinely used in hydrogeochemistry including:

- Ternary diagrams
- Piper plots
- Modified Piper plots
- Durov plots
- Stiff diagrams
- Maucha diagrams
- Multirectangular diagrams
- Schoeller diagrams
- Radial diagrams
- ilr compositional plots

The plots used in the analysis of the chemical composition of geothermal fluids include:

- Custom ternary diagrams
- Giggenbach plot

The plots used in the analysis of the chemical composition of geothermal gases include:

- Custom ternary diagrams
- Grid geothermometers

The package includes functions to calculate several geothermometers including:

- Silica geothermometers
- Fournier-Potter geothermometer
- Na-K geothermometers
- Na-K-Ca geothermometer
- K-Mg geothermometer
- Mg-Li geothermometer

In addition, the package includes a set of functions to study mixing of different types of groundwater using the M3 algorithm.

## Installation 

The GQAnalyzer package is not available on CRAN and therefore it must be installed from github using:

```r
devtools::install_github("khaors/GQAnalyzer")
```

## Usage

The GQAnalyzer package does not require compilations and once installed it can be directly used by loading it: 

```r
library(GQAnalyzer)
```

## Package Documentation

The package documentation can be accesed [here](https://khaors.github.io/GQAnalyzer/) 

## Shiny App

This package includes a shiny app called _GQAnalyzer\_gui_ designed to help in the interpretation of the chemical composition of groundwater, geothermal fluids and/or geothermal gases samples. This app can be called using:

```r
GQAnalyzer_gui()
```

## Contact
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/khaors/GQAnalyzer/issues)

Want to help or have questions? Contact me directly, use an issue, fork me or submit a pull request.

## License

GPL
