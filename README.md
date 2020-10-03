# AEM-ISP
Acousto-Electric-Modulation Inverse Source Problem

## Before you run
Modify the Makefile.in in ``femm`` to get the path to Matlab correct. Then install the packages (Ubuntu for example)

```
apt-get install libtriangle-dev
apt-get install libmetis-dev
```

## How to run
Before the first run in Matlab, change directory to aem-isp.
```
startup; cd femm; make; cd ..
```
After the compilation finishes, then run the demo.
```
demo;
```

## How to use
Just change the functions in ``utils``.