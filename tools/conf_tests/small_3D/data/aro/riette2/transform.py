#!/usr/bin/env python3

import epygram
epygram.init_env()

with epygram.formats.resource('ICMSHFCSTINIT_l15', 'a') as r:
    for fid in [fid for fid in r.listfields() if fid.endswith('HUMI.SPECIFI')]:
        f = r.readfield(fid)
        f.setdata(f.getdata() * 1.5)
        r.writefield(f)

