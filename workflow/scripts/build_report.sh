#!/bin/bash

# copy files
cp $1 packages/visualizer/src/assets/results.json
cp $2 packages/visualizer/src/assets/loss.json
cp $3 packages/visualizer/src/assets/regions.json
cp $4 packages/visualizer/src/assets/amplicons.json

# move to visualizer path
CURRENT_DIR=$(pwd)
cd packages/visualizer

# install dependencies
echo "Installing dependencies"
npm install
npm install --save-dev gulp gulp-inline

# build
echo "Building visualizer"
npm run-script build
npx gulp

# replace media="print" with type="text/css"
# This is a hack to make the css work in the browser
echo "Modifying CSS"
sed -i 's/media="print"/type="text\/css"/g' single-dist/index.html
echo "Done building visualizer"

# copy to output
cd $CURRENT_DIR
mv packages/visualizer/single-dist/index.html $5

# cleanup
echo "Cleaning up"
cd packages/visualizer
rm dist -rf
rm single-dist -rf
