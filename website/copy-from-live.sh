#!/usr/bin/env bash
#
# After making changes to the live website, we can copy them here.

src=villani:/srv/projects/ircolitis

command rsync -vr $src/{index.html,img,css} ./

mkdir -p app

command rsync -vr $src/app/{index.html,css,ext,js,make-symlinks.sh} ./app/

mkdir -p gene-contrasts

command rsync -vr $src/app/{index.html,css,ext,js,make-symlinks.sh} ./gene-contrasts/

