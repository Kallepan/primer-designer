rule install_packages:
    conda: "../envs/nodejs.yml"
    output: "packages/visualizer/node_modules"
    input: "results/{species}.{pool}.set.json"
    shell: "cd packages/visualizer && npm install && npm install --save-dev gulp && gulp gulp-inline"

rule build_report:
    conda: "../envs/nodejs.yml"
    input: "packages/visualizer/node_modules"
    output: "packages/visualizer/single-dist/index.html"
    shell: "cd packages/visualizer && npm run-script build && npx gulp"

rule copy_to_target:
    conda: "../envs/nodejs.yml"
    input: "packages/visualizer/single-dist/index.html"
    output: "results/summary.html"
    shell: 
    """
        cp {input} {output}
    """

TODO