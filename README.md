# Stibium

*Stibium* is a computer-aided design (CAD) tool from a parallel universe
in which CAD software evolved from Lisp machines rather than drafting tables.

It is a continuation of [Antimony](https://github.com/mkeeter/antimony) by
[Matt Keeter](https://mattkeeter.com/projects/antimony) — itself a spiritual
successor to [kokopelli](https://github.com/mkeeter/kokopelli) by way of
[fabserver](http://kokompe.cba.mit.edu). The name keeps the lineage: Sb is
antimony's element symbol, and every `.sb` model file already knew.

Stibium picks up where Antimony's maintenance mode left off. So far:
mesh-optimized exports (error-bounded simplification), chamfer and fillet
CSG operations, feature detection rehabilitated and enabled by default,
and assorted crash fixes. The roadmap lives in [TODO.md](TODO.md) and
[doc/LIBRARY-ROADMAP.md](doc/LIBRARY-ROADMAP.md).

## Try it
To get started, [build from source](BUILDING.md) (Mac and Linux).

There is also a community-supported package for [Fedora 22 or later](https://admin.fedoraproject.org/pkgdb/package/antimony/):
```
dnf install antimony
```
(or `dnf install antimony --enablerepo=updates-testing` to get testing builds)

## Support

If you hit an issue or have a question, [open an issue](https://github.com/nsfm/stibium/issues).

For Antimony history and design rationale, see [Matt Keeter's writeup](https://mattkeeter.com/projects/antimony).

## License

Stibium's new work is released under the [GPLv3](LICENSE).

Code inherited from Antimony remains under its original MIT License —
see [LICENSE.antimony](LICENSE.antimony).

Copyright (c) 2013-2015 Matthew Keeter and other contributors (Antimony)

Antimony includes code from [kokopelli](https://github.com/mkeeter/kokopelli), which is  
© 2012-2013 Massachusetts Institute of Technology  
© 2013 Matthew Keeter
