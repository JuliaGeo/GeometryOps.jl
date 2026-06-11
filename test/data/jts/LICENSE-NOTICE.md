# JTS test data license notice

The XML files in this directory tree (`test/data/jts/general/`,
`test/data/jts/validate/`) are vendored, unmodified, from the
[JTS Topology Suite](https://github.com/locationtech/jts) test resources
(`modules/tests/src/test/resources/testxml/`).

- Upstream repository: https://github.com/locationtech/jts
- Upstream commit: `123a182e6e5a9cc8caed8ff037e4f824a5ce74ee` (2026-03-05)
- Vendored on: 2026-06-11

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
