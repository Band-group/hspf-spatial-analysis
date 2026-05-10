#!/bin/bash
DATE=`date +'%Y-%m-%d'`
DEST=gav@linux.well.ox.ac.uk:/home/gav/public_html/projects/pfsa/simulation/${DATE}

echo "++ Building site.."
npm run build
echo "++ Deploying site to ${DEST}..."
rsync -av dist/ ${DEST}

echo "++ To make this 'live' please update the relevant symlink."

