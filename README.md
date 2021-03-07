# gtfsrouter <a href='https://atfutures.github.io/gtfs-router/'><img src='man/figures/gtfsrouter.png' align="right" height=210 width=182/></a>

[![R build
status](https://github.com/atfutures/gtfs-router/workflows/R-CMD-check/badge.svg)](https://github.com/atfutures/gtfs-router/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ATFutures/gtfs-router/branch/master/graph/badge.svg)](https://codecov.io/gh/ATFutures/gtfs-router)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/gtfsrouter)](https://cran.r-project.org/package=gtfsrouter)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/gtfsrouter?color=orange)](https://cran.r-project.org/package=gtfsrouter)

**R** package for routing with [GTFS (General Transit Feed
Specification)](https://developers.google.com/transit/gtfs/) data. See
[the website](https://atfutures.github.io/gtfs-router/) for full
details.

## Installation

To install:

``` r
remotes::install_github("atfutures/gtfs-router")
```

You can install latest stable version of `gtfsrouter` from CRAN with:

``` r
install.packages("gtfsrouter") # current CRAN version
```

Alternatively, current development versions can be installed using any
of the following options:

``` r
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/gtfsrouter")
remotes::install_bitbucket("atfutures/gtfsrouter")
remotes::install_gitlab("atfutures1/gtfsrouter")
remotes::install_github("ATFutures/gtfsrouter")
```

To load the package and check the version:

``` r
library(gtfsrouter)
packageVersion("gtfsrouter")
```

    ## [1] '0.0.4'

## Main functions

The main functions can be demonstrated with sample data included with
the package from Berlin (the Verkehrverbund Berlin Brandenburg, or VBB).
GTFS data are always stored as `.zip` files, and these sample data can
be written to local storage with the function `berlin_gtfs_to_zip()`.

``` r
berlin_gtfs_to_zip()
tempfiles <- list.files (tempdir (), full.names = TRUE)
filename <- tempfiles [grep ("vbb.zip", tempfiles)]
filename
```

    ## [1] "/tmp/RtmpLlSA7K/vbb.zip"

For normal package use, `filename` will specify the name of the local
GTFS data stored as a single `.zip` file.

### gtfs\_route

Given the name of a GTFS `.zip` file, `filename`, routing is as simple
as the following code:

``` r
gtfs <- extract_gtfs (filename)
gtfs <- gtfs_timetable (gtfs) # A pre-processing step to speed up queries
gtfs_route (gtfs,
            from = "Schonlein",
            to = "Berlin Hauptbahnhof",
            start_time = 12 * 3600 + 120) # 12:02 in seconds
```

| route\_name | trip\_name    | stop\_name                      | arrival\_time | departure\_time |
|:------------|:--------------|:--------------------------------|:--------------|:----------------|
| U8          | S+U Wittenau  | U Schonleinstr. (Berlin)        | 12:09:00      | 12:09:00        |
| U8          | S+U Wittenau  | U Kottbusser Tor (Berlin)       | 12:11:00      | 12:11:00        |
| U8          | S+U Wittenau  | U Moritzplatz (Berlin)          | 12:13:00      | 12:13:00        |
| U8          | S+U Wittenau  | U Heinrich-Heine-Str. (Berlin)  | 12:14:30      | 12:14:30        |
| U8          | S+U Wittenau  | S+U Jannowitzbrucke (Berlin)    | 12:15:30      | 12:15:30        |
| S3          | S Spandau Bhf | S+U Jannowitzbrucke (Berlin)    | 12:20:54      | 12:21:24        |
| S3          | S Spandau Bhf | S+U Alexanderplatz Bhf (Berlin) | 12:22:54      | 12:23:42        |
| S3          | S Spandau Bhf | S Hackescher Markt (Berlin)     | 12:24:54      | 12:25:24        |
| S3          | S Spandau Bhf | S+U Friedrichstr. Bhf (Berlin)  | 12:26:54      | 12:27:42        |
| S3          | S Spandau Bhf | S+U Berlin Hauptbahnhof         | 12:29:36      | 12:30:12        |

### gtfs\_isochrone

*This function will soon be deprecated, to be replaced by
`gtfs_traveltimes()`.*

Isochrones from a nominated station - lines delineating the range
reachable within a given time - can be extracted with the
`gtfs_isochrone()` function, which returns a list of all stations
reachable within the specified time period from the nominated station.

``` r
gtfs <- extract_gtfs (filename)
gtfs <- gtfs_timetable (gtfs) # A pre-processing step to speed up queries
x <- gtfs_isochrone (gtfs,
                     from = "Schonlein",
                     start_time = 12 * 3600 + 120,
                     end_time = 12 * 3600 + 720) # 10 minutes later
```

The function returns an object of class `gtfs_isochrone` containing
[`sf`](https://github.com/r-spatial/sf)-formatted sets of start and end
points, along with all intermediate (“mid”) points, and routes. An
additional item contains the non-convex (alpha) hull enclosing the
routed points. This requires the packages
[`geodist`](https://github.com/hypertidy/geodist),
[`sf`](https://cran.r-project.org/package=sf),
[`alphahull`](https://cran.r-project.org/package=alphahull), and
[`mapview`](https://cran.r-project.org/package=mapview) to be installed.
Isochrone objects have their own plot method:

``` r
plot (x)
```

<img src="man/figures/isochrone.png" width = "80%"/>

The isochrone hull also quantifies its total area and width-to-length
ratio.

## Additional Functionality

There are many ways to construct GTFS feeds. For background information,
see [`gtfs.org`](http://gtfs.org), and particularly their [GTFS
Examples](https://docs.google.com/document/d/16inL5BVcM1aU-_DcFJay_tC6Ni0wPa0nvQEstueG5k4/edit).

Feeds may include a “frequencies.txt” table which defines “service
periods”, and overrides any schedule information during the specified
times. The `gtfsrouter` package includes a function,
[`frequencies_to_stop_times()`](https://atfutures.github.io/gtfs-router/reference/frequencies_to_stop_times.html),
to convert “frequencies.txt” tables to equivalent “stop\_times.txt”
entries, to enable the feed to be used for routine.

Feeds may also omit a “transfers.txt” table which otherwise defines
transfer abilities and times between different services. Feeds without
this table can generally not be used for routing, and they exclude the
possibility of transferring between multiple services. The `gtfsrouter`
package also includes a function,
[`gtfs_transfer_table()`](https://atfutures.github.io/gtfs-router/reference/gtfs_transfer_table.html),
which can calculate a transfer table for a given feed, with transfer
times calculated either using straight-line distances (the default), or
using more realistic times routed through the underlying street network.

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) following the
[all-contributors](https://allcontributors.org) specification.
Contributions of any kind are welcome!

### Code

<table>
<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/AlexandraKapp">
<img src="https://avatars.githubusercontent.com/u/18367515?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/commits?author=AlexandraKapp">AlexandraKapp</a>
</td>
<td align="center">
<a href="https://github.com/stmarcin">
<img src="https://avatars.githubusercontent.com/u/11378350?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/commits?author=stmarcin">stmarcin</a>
</td>
<td align="center">
<a href="https://github.com/polettif">
<img src="https://avatars.githubusercontent.com/u/17431069?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/commits?author=polettif">polettif</a>
</td>
</tr>
</table>

### Issues

<table>
<tr>
<td align="center">
<a href="https://github.com/tbuckl">
<img src="https://avatars.githubusercontent.com/u/98956?u=9580c2ee3c03cbbe44ac8180b0f6a6725b0415f0&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Atbuckl">tbuckl</a>
</td>
<td align="center">
<a href="https://github.com/sridharraman">
<img src="https://avatars.githubusercontent.com/u/570692?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Asridharraman">sridharraman</a>
</td>
<td align="center">
<a href="https://github.com/tuesd4y">
<img src="https://avatars.githubusercontent.com/u/13107179?u=cfcc7852d1bed6e2b17fa3f985cebf743c43b299&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Atuesd4y">tuesd4y</a>
</td>
<td align="center">
<a href="https://github.com/luukvdmeer">
<img src="https://avatars.githubusercontent.com/u/26540305?u=c576e87314499815cbf698b7781ee58fd1d773e2&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Aluukvdmeer">luukvdmeer</a>
</td>
<td align="center">
<a href="https://github.com/Robinlovelace">
<img src="https://avatars.githubusercontent.com/u/1825120?u=461318c239e721dc40668e4b0ce6cc47731328ac&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3ARobinlovelace">Robinlovelace</a>
</td>
<td align="center">
<a href="https://github.com/orlandoandradeb">
<img src="https://avatars.githubusercontent.com/u/48104481?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Aorlandoandradeb">orlandoandradeb</a>
</td>
<td align="center">
<a href="https://github.com/Maxime2506">
<img src="https://avatars.githubusercontent.com/u/54989587?u=6d2c848ee0c7d8a2841d47791c30eec1cab35470&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3AMaxime2506">Maxime2506</a>
</td>
</tr>
<tr>
<td align="center">
<a href="https://github.com/chinhqho">
<img src="https://avatars.githubusercontent.com/u/47441312?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Achinhqho">chinhqho</a>
</td>
<td align="center">
<a href="https://github.com/federicotallis">
<img src="https://avatars.githubusercontent.com/u/25511806?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Afedericotallis">federicotallis</a>
</td>
<td align="center">
<a href="https://github.com/rafapereirabr">
<img src="https://avatars.githubusercontent.com/u/7448421?u=9a760f26e72cd66150784babc5da6862e7775542&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Arafapereirabr">rafapereirabr</a>
</td>
<td align="center">
<a href="https://github.com/loanho23">
<img src="https://avatars.githubusercontent.com/u/48426365?u=36727e1ed27b3b6206fb922c47544ef249fad83d&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Aloanho23">loanho23</a>
</td>
<td align="center">
<a href="https://github.com/dcooley">
<img src="https://avatars.githubusercontent.com/u/8093396?u=2c8d9162f246d90d433034d212b29a19e0f245c1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Adcooley">dcooley</a>
</td>
<td align="center">
<a href="https://github.com/dhersz">
<img src="https://avatars.githubusercontent.com/u/1557047?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ATFutures/gtfs-router/issues?q=is%3Aissue+author%3Adhersz">dhersz</a>
</td>
</tr>
</table>
<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
