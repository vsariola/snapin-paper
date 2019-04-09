# snapin-paper

This repository contains the code related to the paper:

Veikko Sariola "Analytical expressions for spring constants of capillary bridges and snap-in forces of hydrophobic surfaces"

If you find the code useful, consider citing that paper.

## What is it?

The repository contains scripts to recreate the plots in that paper. It uses following submodules:
- kenmotsu-drop-simulator. This module, written by the author, is the simulator for solving the capillary bridge shapes numerically.
- sideaxes-m. This module, written by the author, is used to create nice ticks, labels and scales for the plots.
- export_fig. This submodule is used make high quality pngs from the plots. Written by Oliver J. Woodford & Yair M. Altman.

## How to recreate the figures from the paper

### Prerequisites

- The code was developed on Matlab R2017b, running on Windows. May or may not work on other versions of Matlab.
- Latest version of the code. Can be downloaded from https://github.com/vsariola/drop-simulator
- To export figures into pdfs, export_fig needs Ghostscript.
- Some of the simulations use parallelization to speed them up; if you do not have Parallel Computing Toolbox, replace the parfor in par_fun.m with an ordinary for.


### Figure 1

1. Change Matlab working directory to `Figure Intro\`
2. Run `simulate_constant_angle.m`
3. Run `simulate_constant_radius.m`
4. Run `plot_constant_angle.m`
5. Run `plot_constant_radius.m`

The figures are still missing some labels; the final composition was done in Inkscape (`Figure_Intro.svg`)

### Figure 2

1. Change Matlab working directory to `Figure Schematic\`
2. Run `plot_schematic.m`

The figure was heavily edited in Inkscape; the final version is shown in `Figure_Schematic.svg`.

### Figure 3

1. Change Matlab working directory to `Figure Analytical vs Numerical\`
2. Run `simulate_wa_vs_snapin.m`
3. Run `simulate_r_vs_snapin.m`
4. Run `plot_wa_vs_snapin.m` 

Additionally, you can run `find_limit_of_accuracy.m` to get the numbers reported as the accuracy of the approximation in the paper.

### Figure 4

1. Change Matlab working directory to `Figure Experimental Force Distance`
2. Run `process_experimental_fd.m`
3. Run `simulate_fd.m`
4. Run `plot_experimental_fd.m`

### Figure 5

1. Change Matlab working directory to `Figure Model vs Experimental`
2. Run `simulate_experimental_comparison.m`
3. Run `plot_model_vs_experimental.m`

### Figure 6

1. Change Matlab working directory to `Figure Simulated Experiment`
2. Run `simulate_experiment.m`
3. Run `plot_simulated_experiment.m`

## Licensing

See LICENSE.

The code uses the export_fig library, see export_fig/LICENSE for its license.

kenmotsu-drop-simulator uses cumquad.m, see kenmotsu-drop-simulator/cumquad_license.txt for its license.


