# JTS test data license notice

The XML files in this directory tree (`test/data/jts/general/`,
`test/data/jts/validate/`, `test/data/jts/overlay/`,
`test/data/jts/overlay_robust/`) are vendored, unmodified, from the
[JTS Topology Suite](https://github.com/locationtech/jts) test resources
(`modules/tests/src/test/resources/testxml/`).

- Upstream repository: https://github.com/locationtech/jts
- Upstream commit: `123a182e6e5a9cc8caed8ff037e4f824a5ce74ee` (2026-03-05)
- Vendored on: 2026-06-11 (relate suites), 2026-07-28 (overlay suites)

JTS is dual-licensed under the Eclipse Public License 2.0 (EPL 2.0) and the
Eclipse Distribution License 1.0 (EDL 1.0, a BSD-style license). These files
are redistributed here under those terms. See:

- https://github.com/locationtech/jts/blob/master/LICENSE_EPLv2.txt
- https://github.com/locationtech/jts/blob/master/LICENSE_EDLv1.txt

File provenance within the upstream `testxml/` directory:

- `general/TestRelate{PP,PL,PA,LL,LA,AA}.xml`, `general/TestBoundary.xml` → `test/data/jts/general/`
- `misc/TestRelateEmpty.xml`, `misc/TestRelateGC.xml` → `test/data/jts/general/`
- `validate/TestRelate*.xml` → `test/data/jts/validate/`
- `robust/TestRobustRelate.xml`, `robust/TestRobustRelateFloat.xml` → `test/data/jts/validate/`
- `general/TestOverlay{AA,LA,LL,PA,PL,PP,Empty}.xml`,
  `general/TestNGOverlay{A,L,P,Empty,GC}.xml` → `test/data/jts/overlay/`
- `misc/TestOverlay.xml` → `test/data/jts/overlay/TestOverlayMisc.xml`
  (renamed only to avoid a basename collision with `general/TestOverlay*.xml`;
  contents unmodified)
- `robust/overlay/*.xml` up to 140 kB → `test/data/jts/overlay_robust/` (44 of
  the 51 files, 972 kB). The seven omitted files are the multi-megabyte ones —
  `TestOverlay-stmlf.xml` (6.4 MB), `TestOverlay-geos-997-union-slow.xml`
  (4.5 MB), `TestOverlay-geos-1034.xml`, `TestOverlay-pg-4738.xml`,
  `TestOverlay-geos-600-lines.xml`, `TestOverlay-isochrone.xml`,
  `TestOverlay-geos-392-lines.xml` — which are too heavy to carry in-tree; run
  them out of a JTS checkout if a deeper sweep is ever wanted.

Not vendored, by design: the `*Prec.xml` overlay suites
(`general/TestOverlay{AA,LA,LL,PL}Prec.xml`,
`general/TestNGOverlay{A,L,P}Prec.xml`). They declare a fixed
`<precisionModel scale=...>`, and fixed-precision / snap-rounding overlay is
permanently out of scope for this engine (OverlayNG port design §0).
