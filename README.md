# optimal sensing

This is the code for the paper `Data-Driven Sensor Placement for Fluid Flows`.
The description of the folders is as follows.

Folder       | DESCRIPTION
-----------------|-------------
`src`        | contains codes to perform model learning and sensor parameter optimization
`ginzburg_landau` | codes to generate the state transition matrix for the linearized complex Ginzbug-Landau system
`scripts`         | Matlab scripts to generate figures in the paper
`figures`          | Figures in the paper
`fortran_wrapper`  | Fortran code that performs observer-based feedback flow control on the flow over an inclined flat plate

To make all the matlab scripts accessible, run `import_all.m`.
The flow over an inclined flat plate is simulated using the following code -- https://github.com/cwrowley/ibpm.
