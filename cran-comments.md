## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new release. It has been checked on a local Ubuntu install with R 4.4.1, win-builder release and dev, and osx using Rhub. 

* checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      libs   5.2Mb

* checking compiled code ... NOTE
  File ‘xactonomial/libs/xactonomial.so’:
    Found non-API calls to R: ‘BODY’, ‘CLOENV’, ‘DATAPTR’, ‘ENCLOS’,
      ‘FORMALS’

The size is due to the vendoring of Rust libraries. 

I'm not sure if the last note is a deal-breaker. It comes from a dependency on libR-sys. My understanding is that the C API entry points are still subject to change, and the libR-sys developers are actively monitoring the situation and will update as needed to be compliant. 

