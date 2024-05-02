#!/bin/bash -ex
#/hlilab/hli/gafcall/chm13.reg.bed 


gafcall.js join graph-sv hg38-sv | gafcall.js chm13-sv - > filtered-sv

../gafcall/js/gafcall.js view -C -b chm13.reg.bed /hlilab/hli/gafcall/pair/chm13/COLO829T.hifi1.pair.gsv
#../gafcall/js/gafcall.js view -C -b chm13.reg.bed /hlilab/hli/gafcall/pair/chm13/COLO829T.hifi1.only.gsv

#gafcall.js join graph-sv hg38-sv | gafcall.js chm13-sv - > filtered-sv
