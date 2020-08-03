# TPG4155_Project2018
The project title is "Pragmatic Modeling of a Full Waveform Inversion using Finite Difference Method to Solve Wave Equation". You can see the full project report [here](https://github.com/fdlberylian/TPG4155_Project2018/blob/master/TPG4155_Fall_Project_2018__10003%2610010.pdf).

List of Code Scripts
1) calculateGradient : Function to calculate gradient between test model and true model.
2) calculateStepLength : Function to calculate step length used to minimize error, using "blind" steepest descent as the numerical optimization method.
3) calculateStepLength_v2 : Function to calculate step length used to minimize error, using golden search as the numerical optimization method.
4) rickerWave : Function to give source signature, based on time derivative of Ricker wavelet.
5) runFWI : This is the main script to run the Full Waveform Inversion program.
6) solveWaveEqn : Function to solve acoustic wave equation numerically using finite difference method.
7) taperGradient : Function to taper the gradient near receivers using sine tapering.
