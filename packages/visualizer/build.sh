#!/bin/bash
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
