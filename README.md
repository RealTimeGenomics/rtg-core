# RTG Core Non-commercial

Copyright (c) 2014 Real Time Genomics Ltd

This software is provided under the Real Time Genomics Ltd Software
Licence Agreement for Academic Non-commercial Research Purposes
only. See [LICENSE.txt](LICENSE.txt)

If you require a license for commercial use or wish to purchase
commercial support contact us via info@realtimegenomics.com.

---

## Prerequisites

* Java 1.7 or later
* apache ant 1.9 or later

## Check out source code for both RTG Tools and RTG Core

Building RTG Core requires the source code to both RTG Tools and RTG
Core, so you must clone both repositories:

    $ git clone https://github.com/RealTimeGenomics/rtg-tools.git
    $ git clone https://github.com/RealTimeGenomics/rtg-core.git
    $ cd rtg-core

To update, you will need to perform a `git pull` on both
repositories. Advanced users may use `git subtree` to have both RTG
Tools and RTG Core within a single repository (it's what we use during
development).

## Compile / run unit tests

    $ ant runalltests

## Building RTG Core

To build the RTG Core Non-Commercial package which can be locally
installed and run:

    $ ant zip-nojre

This will create an installation zip file under `dist`.

## Installation

Uncompress the installation zip:

    $ cd /where/you/want/to/install/
    $ unzip /path/to/rtg-core/dist/rtg-core-VERSION-nojre.zip

Follow the instructions contained in the `README.txt`. This build will
use the system Java by default, so ensure it is Java 1.7 or later.

You can use the commands in RTG Core to format your own reference
datasets, or download common
[pre-formatted references](http://realtimegenomics.com/news/pre-formatted-reference-datasets/)
from our website.

For a nice demonstration of the features of RTG Core for sex and
pedigree aware mapping and variant calling on data generated from
scratch using the RTG Core simulation tools, run the `demo-family.sh`
script contained in the scripts of the installation directory:

    $ ./scripts/demo-family.sh $PWD/rtg

A PDF user manual is included in the installation directory or can be
[viewed online](installer/resources/core/RTGOperationsManual.pdf).

## Release history

See [doc/ReleaseNotes.txt](doc/ReleaseNotes.txt) for full release
history details.

## Support

An
[rtg-users](https://groups.google.com/a/realtimegenomics.com/forum/#!forum/rtg-users)
discussion group is now available for general questions, tips, and
other discussions.

To be informed of new software releases, subscribe to the low-traffic
[rtg-announce](https://groups.google.com/a/realtimegenomics.com/forum/#!forum/rtg-announce)
group.

Some email support will be available via `support@realtimegenomics.com`
which we will do our best to fulfill. If you have specific support or
service requirements, talk to us about our support offerings which
include guaranteed response times, training, project advice, priority
bug fixes and feature requests.

