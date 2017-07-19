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

## INSTALL

Please follow this simple instructions.

### Requirements

You just need a working [ROOT](http://root.cern.ch/) installation (ROOT v 5.34.36 or higher is the recommended version). Then you just execute <kbd>make</kbd> at the main directory:

```bash
make
```

It will create the <kbd>${LAGO_EASY}</kbd> environment variable and modify the <kbd>$HOME/.bashrc</kbd> file. Please source it again before to continue.

Please report any installation problems please [email us](mailto:lago@lagoproject.org).

### Usage

#### Modifying geometry

The detector geometry is defined in the <kbd>src/Constants.h</kbd> file. A symbolink link named <kbd>configs.h</kbd> is provided. To modify just edit this file:

```bash
vim configs.h
```

and modify it by changing the values of STATION_RADIUS, STATION_HEIGHT and the number and position of PMTs (NPM X_PM and Y_PM). For example, for one of our detectors deployed at the LAGO Bariloche site, the configuring block is:

``` cpp
//DETECTOR CONSTANTS
const double STATION_RADIUS     = 0.78; //(m)
const double STATION_HEIGHT     = 1.54; //(m)
const int  NPM                  = 1;
const double X_PM[NPM]          = {0.}; // (0,0) is the center of the tank roof
const double Y_PM[NPM]          = {0.}; // (0,0) is the center of the tank roof
const double RAD_PM             = .1477;//(m)
const double HC_PM              = .0776;//(m)
const double TOP_FACT           = 1.;  //white top
```

If you change detector geometry, you will need to recompile the whole package. Just execute <kbd>make</kbd> at the parent directory:

```bash
 cd ${LAGO_EASY}/src
 make
```

#### Running

The simulation execution is governed by the <kbd>default.inp</kbd> file. Just edit if following the examples given in that file. It is straightforward. Please be sure to be in CALIB mode (SHOWER mode is deactivated in LAGO code). Then, for run the simulation just use

```bash
cd ${LAGO_EASY}
make run
```

in the parent directory. It will provide a <kbd>.root</kbd> file with a specific name depending on the parameters defined in the <kbd>default.inp</kbd> file. It is then converted into ASCII, extracting the FADC traces, producing a <kbd>.dat</kbd> ASCII file.

#### Data analysis

The output data is quite simple. In the resulting ASCII file, the expected ADC values for the first 100 temporal bins of the FADC traces (~2.5 us) are given for each injected particle. A simple awk script to produce a charge histogram could be: 
```bash
for i in *.dat; do 
	u=$(basename $i .dat)
	echo $i
	awk '{t[$1]+=$4}END{for (i in t) {print t[i] }}'  $i | sort -n | uniq -c | awk '{if ($2) {print $2,$1}}' > $u.hst
```
producing a file containing a two column file: <kbd>charge counts</kbd>.

#### Housekeeping

We include three leves of data and binary cleaning: 
```bash
make clean		# clean binaries and libraries
make simclean   # erase simulated data from previous runs
make fullclean  # executes clean and simclean directives 
```

### Additional package

An additional tool, <kbd>Trace</kbd> is included. It can be used to calculate the trace lenght distribution of muons going thought the detector. To configure it, two blocks inside the main code should be modified. So: 
```bash
vim src/trace.cc
```
and then, modify these two blocks according to your needs:
```cpp
// config
const double deg2rad  = M_PI / 180.;
const double rad2deg  = 180. / M_PI;
const bool   internal = true; // use internal generator (true) or read muons from external file
const bool   extSeeds = true;  // use rand() to generate the seeds or just set them as -1
const double muons    =  1e6;  // number of showers to generate if internal is enabled
const double units    =  10.;  // in milimeters
```
and
```cpp
// guane-3
const double tankR = 515 / units; // detector radius
const double tankH = 500. / units; // detector height
#define MU "guane.muo" // muon crossing the detector distribution file
#define TR "guane.trc" // traces distribution file
#define DI "guane.dst" // all muons speed distribution file
```
Then, just compile and run Trace:
```bash
cd $LAGO_EASY
make
./bin/Trace
```
