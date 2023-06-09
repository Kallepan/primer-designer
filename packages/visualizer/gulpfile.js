// Create gulpfile.js
const gulp = require("gulp");
const inline = require("gulp-inline");

gulp.task("default", () => {
  return gulp
    .src("./dist/index.html")
    .pipe(inline())
    .pipe(gulp.dest("./single-dist"));
});

// Install dependencies
// npm install --save-dev gulp gulp-inline
// Build Angular app with ng build
// after that run npx gulp to compile the build output to a single file
// or add this to package.json "build-single": "npm run build && npx gulp"