# NonintegerCollatz

This code is part of a [project](https://blog.doingsciencetostuff.com/2026/07/13/extending-collatz-to-non-integers/) exploring Collatz for non-integers.

## Getting Started

### Dependancies

This code relies on the GNU Multiple Precision Arithmetic Library for
calculation and GNU Plot for making graphs.
Packages for [libgmp](https://gmplib.org) and [gnuplot](http://www.gnuplot.info) should exist for most Linux distributions and
are available via [macports](https://www.macports.org) on macOS.

In [Debian](https://www.debian.org/) Linux:
```text
$ sudo apt-get install libgmp10-dev gnuplot
```

In macOS using [macports](https://www.macports.org):
```text
$ sudo port install gmp gnuplot
```

### Running

```
$ git clone https://github.com/fra99le/NonintegerCollatz.git
$ cd NonintegerCollatz
$ c++ -I/opt/local/include -L/opt/local/lib collatz_w.cpp -lgmp -lgmpxx -o collatz_w 
$ ./collatz_w 33
$ ./plot_w.gplot
```
