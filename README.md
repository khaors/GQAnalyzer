# __GQAnalyzer__: Package to create hydrogeochemical plots
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rves)](https://cran.r-project.org/package=GQAnalyzer)
[![HitCount](http://hits.dwyl.io/khaors/rves.svg)](http://hits.dwyl.io/khaors/GQAnalyzer)

The goal of GQAnalyzer is to create several plots routenely used in hydrogeochemistry including:

- Ternary diagrams
- Piper plots
- Durov plots
- Stiff diagrams
- Maucha diagrams
- Multirectangular diagrams
- Schoeller diagrams
- Radial diagrams

## Installation 

The GQAnalyzer package is not available on CRAN and therefore it must be installed from github using:

```r
devtools::install_github("khaors/GQAnalyzer")
```

## Usage

The pumpingtest package does not require compilations and once installed it can be directly used by loading it: 

```r
library(GQAnalyzer)
```

## Package Documentation

The package documentation can be accesed [here](https://khaors.github.io/pumpingtest/) 

## Shiny App

This package includes a shiny app called _GQAnalyzer\_gui_ designed to help in the interpretation of the Vertical Electric Soundings. This app can be called using:

```r
GQAnalyzer_gui()
```

## Contact
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/khaors/GQAnalyzer/issues)

Want to help or have questions? Contact me directly, use an issue, fork me or submit a pull request.

## License

GPL
