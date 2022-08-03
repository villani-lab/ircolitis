#!/usr/bin/env bash
#
# After making changes to the live website, we can copy them here.

src=/srv/projects/ircolitis

command cp -rf $src/{index.html,tables.html,img,css} ./

mkdir -p app

command cp -rf $src/app/{index.html,css,ext,js,make-symlinks.sh} ./app/

mkdir -p de/data
command cp -rf $src/de/{bundle.js,bundle.js.map,bundle.worker.js,bundle.worker.js.map,index.html,bundle.js.LICENSE.txt,sql-wasm-3896522bd16c305a8c29.wasm,7edd8e1eedab0513cc67d0914307919e.svg} ./de/

