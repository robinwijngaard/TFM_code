#!/bin/bash

[ -w ".Rprofile" ] && rm .Rprofile

/home/robin/Downloads/R-3.1.2/bin/Rscript sessionInfo.R --bootstrap-packrat > setup.log 2>&1

cp packrat/packrat_source/.Rprofile ./
