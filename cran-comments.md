## R CMD check results

0 errors | 0 warnings | 2 notes

* This is an update to a new release submission. It has been checked on a local Ubuntu install with R 4.4.1, win-builder release and dev, and osx using Rhub. I removed the dependence on libR-sys and seems to fix the API entry point warnings. I also fixed the broken url in description. 

* New submission

* checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      libs   5.2Mb

The size is due to the vendoring of Rust libraries. 

