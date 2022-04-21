# Octave/MATLAB implementation of GeographicLib

Contents:
* [Introduction](#introduction)
* [Installation](#installation)
* [Function summary](#function-summary)
* [Changes](#changes)
* [Other links](#other-links)

## Introduction

This toolbox provides native MATLAB implementations of a subset of the
C++ library, GeographicLib.  Key components of this toolbox are
  * Geodesics, direct, inverse, area calculations.
  * Projections, transverse Mercator, polar stereographic, etc.
  * Grid systems, UTM, UPS, MGRS.
  * Geoid lookup, egm84, egm96, egm2008 geoids supported.
  * Geometric transformations, geocentric, local cartesian.
  * Great ellipse, direct, inverse, area calculations.

All the functions are vectorized and so offer speeds comparable to
compiled C++ code when operating on arrays.

Some common features of these functions:
  * Angles (latitude, longitude, azimuth, meridian convergence) are
    measured in degrees.
  * Distances are measured in meters, areas in meters^2.
  * Latitudes must lie in `[-90, 90]`.  Latitudes outside this range
    are treated in NaNs
  * The ellipsoid is specified as `[a, e]`, where `a` = equatorial radius
    and `e` = eccentricity.  The eccentricity can be pure imaginary to
    denote a prolate ellipsoid.
  * Keep `|e| < 0.2` (i.e., `|f| <= 1/50`) for full double precision
    accuracy.

There is some overlap between this toolbox and MATLAB's Mapping
Toolbox.  However, this toolbox offers:
  * better accuracy;
  * treatment of oblate and prolate ellipsoid;
  * guaranteed convergence for geoddistance;
  * calculation of area and differential properties of geodesics;
  * ellipsoidal versions of the equidistant azimuthal and gnomonic
    projections.

## Installation

The Octave/MATLAB packages are available on
[SourceForge](
https://sourceforge.net/projects/geographiclib/files/distrib-Octave)
and
[MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/50605).

Tha Octave package can be installed with the command
```octave
pkg install geographiclib-octave-M.N.tar.gz
pkg load geographiclib
```
Other useful Octave package commands
```octave
pkg list                            % list all packages
ver geographiclib                   % list version
pkg describe -verbose geographiclib % list functions
news geographiclib                  % change log
pkg test geographiclib              % run the tests
pkg unload geographiclib            % remove from path
pkg uninstall geographiclib         % uninstall
```

The MATLAB toolbox can be installed with the command
```octave
matlab.addons.install geographiclib_toolbox-M.N.mltbx
```
Other useful MATLAB toolboxes commands
```octave
t = matlab.addons.installedAddons        % list all toolboxes
n = 1;                                   % the row with geographiclib
matlab.addons.uninstall(t.Identifier(n)) % uninstall
```

## Function summary

### Geodesics
  * [`geoddistance`](inst/gedistance.m) -
    Distance between points on an ellipsoid
  * [`geodreckon`](inst/geodreckon.m) -
    Point at specified azimuth, range on an ellipsoid
  * [`geodarea`](inst/geodarea.m) -
    Surface area of polygon on an ellipsoid

### Projections
  * [`tranmerc_fwd`](inst/tranmerc_fwd.m) -
    Forward transverse Mercator projection
  * [`tranmerc_inv`](inst/tranmerc_inv.m) -
    Inverse transverse Mercator projection
  * [`polarst_fwd`](inst/polarst_fwd.m) -
    Forward polar stereographic projection
  * [`polarst_inv`](inst/polarst_inv.m) -
    Inverse polar stereographic projection
  * [`eqdazim_fwd`](inst/eqdazim_fwd.m) -
    Forward azimuthal equidistant projection
  * [`eqdazim_inv`](inst/eqdazim_inv.m) -
    Inverse azimuthal equidistant projection
  * [`cassini_fwd`](inst/cassini_fwd.m) -
    Forward Cassini-Soldner projection
  * [`cassini_inv`](inst/cassini_inv.m) -
    Inverse Cassini-Soldner projection
  * [`gnomonic_fwd`](inst/gnomonic_fwd.m) -
    Forward ellipsoidal gnomonic projection
  * [`gnomonic_inv`](inst/gnomonic_inv.m) -
    Inverse ellipsoidal gnomonic projection

### Grid systems
  * [`utmups_fwd`](inst/utmups_fwd.m) -
    Convert to UTM/UPS system
  * [`utmups_inv`](inst/utmups_inv.m) -
    Convert from UTM/UPS system
  * [`mgrs_fwd`](inst/mgrs_fwd.m) -
    Convert UTM/UPS coordinates to MGRS
  * [`mgrs_inv`](inst/mgrs_inv.m) -
    Convert MGRS to UTM/UPS coordinates

### Geoid lookup
  * [`geoid_height`](inst/geoid_height.m) -
    Compute the height of the geoid above the ellipsoid
  * [`geoid_load`](inst/geoid_load.m) -
    Load a geoid model

### Geometric transformations
  * [`geocent_fwd`](inst/geocent_fwd.m) -
    Conversion from geographic to geocentric coordinates
  * [`geocent_inv`](inst/geocent_inv.m) -
    Conversion from geocentric to geographic coordinates
  * [`loccart_fwd`](inst/loccart_fwd.m) -
    Convert geographic to local cartesian coordinates
  * [`loccart_inv`](inst/loccart_inv.m) -
    Convert local cartesian to geographic coordinates

### Great ellipses
  * [`gedistance`](inst/gedistance.m) -
    Great ellipse distance on an ellipsoid
  * [`gereckon`](inst/gereckon.m) -
    Point along great ellipse at given azimuth and range

### Utility
  * [`defaultellipsoid`](inst/defaultellipsoid.m) -
    Return the WGS84 ellipsoid
  * [`ecc2flat`](inst/ecc2flat.m) -
    Convert eccentricity to flattening
  * [`flat2ecc`](inst/flat2ecc.m) -
    Convert flattening to eccentricity
  * [`geographiclib_test`](inst/geographiclib_test.m) -
    The test suite for the geographiclib package
  * [`geographiclib_signtest`](inst/geographiclib_signtest.m) -
    Another test suite

### Documentation
  * [`geoddoc`](inst/geoddoc.m) -
    Geodesics on an ellipsoid of revolution
  * [`projdoc`](inst/projdoc.m) -
    Projections for an ellipsoid
  * [`gedoc`](inst/gedoc.m) -
    Great ellipses on an ellipsoid of revolution

## Changes

See the [change log](NEWS).  Releases are tagged in git as, e.g.,
[`v1.52`](../../tree/v1.52), [`v2.0`](../../tree/v2.0), etc.

## Other links
  * [GeographicLib](https://geographiclib.sourceforge.io).
  * Implementations in [other languages](
    https://geographiclib.sourceforge.io/doc/library.html#languages).
  * How to install the [geoid datasets](
    https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst).
  * C. F. F. Karney,
    [Transverse Mercator with an accuracy of a few nanometers](
    https://doi.org/10.1007/s00190-011-0445-3),
    J. Geodesy **85**(8), 475–485 (2011);
    [preprint](https://arxiv.org/abs/1002.1417);
    [addenda](https://geographiclib.sourceforge.io/tm-addenda.html).
  * C. F. F. Karney,
    [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z),
    J. Geodesy **87**(1), 43–55 (2013);
    [addenda](https://geographiclib.sourceforge.io/geod-addenda.html).
