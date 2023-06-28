rule install_packages:
    conda: "../envs/nodejs.yml"
    output: "packages/visualizer/node_modules"
    input: "results/{species}.{pool}.set.json"
    shell: "cd packages/visualizer && npm install && npm install --save-dev gulp && gulp gulp-inline"

rule build_report:
    conda: "../envs/nodejs.yml"
    output: "results/report/summary.html"
    input: "packages/visualizer/node_modules"
    shell: "cd packages/visualizer && npm run-script build && npx gulp"

rule copy_to_target:
    conda: "../envs/nodejs.yml"
    output: "results/report/summary.html"
    input: "results/report/summary.html"

