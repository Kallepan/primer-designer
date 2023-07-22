#!/bin/bash
cd $4

# copy files
cp $1 src/assets/results.json
cp $2 src/assets/loss.json
cp $3 src/assets/regions.json
# cp $4 src/assets/amplicons.json

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
mv single-dist/index.html $5

# cleanup
echo "Cleaning up"

rm dist -rf
rm single-dist -rf
