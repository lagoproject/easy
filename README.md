# ![The LAGO Project](http://lagoproject.org/images/lago-logo-90.png "The LAGO Project") The LAGO EASY suite

| CODENAME			| EASY  |
|-------------------|:------|
| COPYRIGHT			| (C) 2012-Today, The LAGO Project, [lagoproject.org](http://lagoproject.org)|
| LICENSE			| BSD-3-Clause |
| REPOSITORY		| https://github.com/lagoproject |
| CONTACT			| [lago@lagoproject.org](mailto:lago@lagoproject.org)|
| DESCRIPTION		| The LAGO project deploys and operates a network of water Cherenkov detectors (WCD) with a unified data structure. Here you will find the latest stable (prod) version of the detector simulation tools of the LAGO project |
| CONTRIBUTORS		| If you want to contribute to this project please [send us an email](mailto:lago@lagoproject.org)|


| CODE GUIDELINES	|		|
|-------------------|:------|
| FILE ENCODING		| UTF8 (please use <kbd>iconv -f YOUR_ENCODING -t UTF-8 file_to_convert > converted_file </kbd> before to push) |
| LANGUAGE			| English (preferred) |
| INDENT STYLE		| [Stroustrup](http://en.wikipedia.org/wiki/Indent_style#Variant:_Stroustrup) using 1 tab for 4 columns wide. check [here for vim setup](http://tedlogan.com/techblog3.html) |
|					| If you prefer, please use: <kbd>astyle -t4 -A4 -y file_to_convert</kbd> before to push
| VERSIONING		| Sequence-based identifiers, v<version>r<release>. First public release: **v1r0**
| INSTALL			| After installing dependences (see INSTALL), just *make*
| USAGE				| Please visit our [wikipage](http://wiki.lagoproject.org) (internal use only)|

The [Latin American Giant Observatory (LAGO)](http://lagoproject.org) is an extended Astroparticle Observatory at global scale. It is mainly oriented to basic research on three branches of Astroparticle physics: the Extreme Universe, Space Weather phenomena, and Atmospheric Radiation at ground level.

# INSTALL

Please follow this simple instructions.

## Requirements

You just need a working [ROOT](http://root.cern.ch/) installation (ROOT v 5.34.36 or higher is the recommended version). Then you just 
 make 
from this directory. It will create the
 LAGO_EASY
environment variable and modify the *.bashrc* file. Please source it again before to continue. 

Please report any installation problems please [email us](mailto:lago@lagoproject.org).

## Usage

### Modifying geometry

The detector geometry is defined in the *Constants.h* file, located at the *src* directory. To modify it:
 cd $LAGO_EASY
 vim src/Constants.h

and modify it by changing the values of STATION_RADIUS, STATION_HEIGHT and the number and position of PMTs (NPM X_PM and Y_PM). For example, for one of our detectors deployed at the LAGO Bariloche site, the configuring block is:
``` cpp
const double STATION_RADIUS     = 0.70; //(m)
const double STATION_HEIGHT     = 1.46; //(m)
const int  NPM                  = 1;
const double X_PM[NPM]          = {0.}; // (0,0) is the center of the tank roof
const double Y_PM[NPM]          = {0.}; // (0,0) is the center of the tank roof
```

Then recompile it
 cd $LAGO_EASY/src
 make

## Running

Go to the data directory,

 cd $LAGO_EASY/data

and edit the *default.inp* file. It is straightforward. Be sure to be in CALIB mode (SHOWER mode is deactivated in LAGO code). Then, just run it from the data directory using *make*:

 make

It will provide a **root** file with a specific name depending on the parameters defined in the *default.inp* file. It is then converted in ASCII, extracting the FADC traces, producing a *.dat* ASCII file.
~                                                                          
