the ilupack code is a mess!!

compiling it yourself seems very involved...

## make process
we do the following now
1. get ilupack archive files from somewhere
  - either there is code in lib/ilupack (you can download it manually from http://ilupack.tu-bs.de/), then we use those archives
  [ NOT IMPLEMENTED - or use the ones from /home/bolle/ilupackV2.2/lib/altix64]
2. combine the archives that we actually need
  - extract: libamd.a libSilupack.a  libZilupack.a libDilupack.a libmetis.a libsparspak.a
  - combine them into libilupack.a
  [- remove extracted object files]


## Versions
Mr Jungemann said that mr. Bollhoefer said: V2.2 and V2.4 dont really differ.
just some new matlab interfaces..

thus, we use the V2.2
