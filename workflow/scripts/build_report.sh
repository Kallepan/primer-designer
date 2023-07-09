#!/bin/bash
cd $3

# copy files
cp ../../$1 src/assets/results.json
cp ../../$2 src/assets/regions.json

# install dependencies
echo "Installing dependencies"
npm install
npm install --save-dev gulp gulp-inline

# cleanup
echo "Cleaning up"
rm dist -rf
rm single-dist -rf

# build
echo "Building visualizer"
npm run-script build
npx gulp

# replace media="print" with type="text/css"
# This is a hack to make the css work in the browser
echo "Modifying CSS"
sed -i 's/media="print"/type="text\/css"/g' single-dist/index.html
echo "Done building visualizer"

mv single-dist/index.html ../../$4