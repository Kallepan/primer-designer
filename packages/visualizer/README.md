# Visualizer

## Description

Visualizer is a TypeScript project built on the Angular 16 framework that visually represents the outcomes of a primer pipeline. The pipeline generates primers for various regions of an organism, producing a bundled HTML report. Additionally, a bed file output containing primers and amplicons separate from visualizer is provided, enabling users to load the data into other genome viewers if desired.

## Installation

- Install [Node.js](https://nodejs.org/en/download/)
- Install [Angular CLI](https://cli.angular.io/)
- Put the regions.json and primers.json files in the assets folder
- Run `npm install && npm install gulp gulp-inline` in the visualizer directory
- Run `bash build.sh` in the visualizer directory

## Usage

After installation and building, a folder called single-dist will appear. Open the contained index.html file in a webbrowser. All data from the pipeline is contained in the index.html file, so no server is required to run the visualizer.