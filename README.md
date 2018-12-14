# Numerical N-Link Pendulum Simulator

These files allow the user to create and analyzed simulations of a generalized n-link pendulum or n-bar linkage. The n-link pendulum can be solved using one of three methods: Newton-Euler, Lagrange, or Differential Algebraic Equations. The n-bar linkage can only be solved using Differential Algebraic Equations.

For descriptions symbolic solutions, see comments in `npend_deriver_NE.m`, `npend_deriver_L.m`, `npend_deriver_DAE.m`, and `nlink_deriver_DAE.m`.

- Additional pendulum functionalities (all prompted in `npend.m`):

    -   Standing wave modes can be animated for any method. Works best for many links with small      oscillation amplitudes relative to total length.

    -   Normal modes of oscillation can be found and animated for n links using the Lagrange          method.

    -   Newton-Euler:   Can add friction as a non-conservative damping force.

    -   DAE Pendulum:   Can oscillate base at any frequency, amplitude, and angle.
    
    -   DAE Linkage:    Can oscillate end at any frequency, amplitude, and angle.
                   :    Can add non-conservative restoring force to first link to animate              cyclical motion.            

- Additional analysis functionalities:

    -   Can compare 2 or 3 simulations at a time (`comp.m`).

    -   Can time run n simulations m times to analyze performance (`teval.m`).

## Setup

None needed, repository has all necessary directories and files to begin.

If a simulation is run for a given number of links using a given solution method, the symbolically derived solution files will be uniquely named and saved in a directory corresponding to the solution method. When the same number of links and the same method are specified in any subsequent simulation, the previously derived solution files will be used.

Saved simulations will be sent to the directory `animations/`.

## Running a simulation

Run `npend.m`, or enter `npend` into the `MATLAB` command line.

All of the simulation parameters and initial conditions will be entered as prompted in the command line. For yes or no questions, enter 'Y' or 'y' for yes or 'N' or 'n' for no. Any letter options can be either uppercase or lowercase.

When a prompt is given with options in parentheses, either the value for the prompt can be entered or the letter for an option.

--> EXAMPLE 1: Playing an animation

    >> npend
    > Load animation? y
     > File Name or [L]ist? l
        DAEp100   Lagrange_10_init_90   NE_10_fric10
     > File Name or [L]ist?  NE_10_fric10                       *No directory or extension needed*
    > [R]eplay or [O]ther or [Q]uit? q

--> EXAMPLE 2: Running a simulation

    >> npend
    > Load animation? n
    > Number of Links: 7
    > [N]ewton-Euler, [L]agrange, [D]AE, or [M]odes? N
    [Deriving solutions]
    Elapsed time is 26.061665 seconds.  
    > Friction ([N]o or value)? 10
    > Link Length or [R]andom or [S]pecified? 3       *Specifying one value sets it for all links*
    > Uniform Initial Angles or [R]andom or [S]pecified? s
     > Initial Angle for Link 1: 90
     > Initial Angle for Link 2: 45
     > Initial Angle for Link 3: 180
     > Initial Angle for Link 4: 0
     > Initial Angle for Link 5: 0
     > Initial Angle for Link 6: 90
     > Initial Angle for Link 7: 180
    > Uniform Initial Angular Velocities or [R]andom or [S]pecified? 0
    > Mass of Links or [P]roportional or [R]andom or [S]pecified? r
     > Link Mass Mean: 1
     > Link Mass Standard Deviation: 0.5
    > Fractional Distance to CoM or [R]andom or [S]pecified? 0.6
    > Acceleration due to Gravity: 1
    > Runtime (or [t]ime/tolerance options): t
     > Runtime: 20
     > Number of timesteps: 400                                      *Default is 10 times runtime*
     > Absolute/Relative Tolerance: 1e-10                                        *Default is 1e-8*
    [Animation]
    [Total Energy Plot]
    > Replay? n
    > Save animation? y
    > File Name: example                                   *'./animations/example.mat' now exists*
    [Writing to file]
    > Run again? n

Other features:

- [M]ode can be selected, this will simulate normal mode oscillation using the Lagrange method via `nmode.m` and `nmode_eig.m`.

- In `npend_DAE.m` and `nlink_DAE.m` (the right-hand side functions for the DAE methods), the user has the option to specify the duration of oscillation, e.g. 1/2 period, 1 period, etc.

- Simulating a pendulum using DAEs in `npend.m`, there is the option to set the driving frequency to activate the standing wave modes of the pendulum. This is intended for the setup in which all links are identical and start at the downward vertical with no initial motion. Additionally, the amplitude of oscillation should be small relative to the total length of the pendulum.

## Additional

- The function `animate.m` is called to animate a simulation from `npend.m`

- The function `npend_deriver_NE.m` is called by `npend.m` to create a solution for an n-link pendulum using the Newton-Euler method.

- The function `npend_deriver_L.m` is called by `npend.m` to create a solution for an n-link pendulum using the Lagrange method.

- The function `npend_deriver_DAE.m` is called by `npend.m` to create a solution for an n-link pendulum using the DAE method.

- The function `nlink_deriver_DAE.m` is called by `npend.m` to create a solution for an n-bar linkage using the DAE method.

- The function `npend_NE.m` (within `Lpend/`) is called by `npend.m` and calls the functions created by `npend_deriver_NE.m` to solve a simulation.

- The function `npend_Lagrange.m` (within `Lpend/`) is called `by npend.m` and calls the functions created by `npend_deriver_L.m` to solve a simulation.

- The function `npend_DAE.m` (within `Lpend/`) is called by `npend.m` and calls the functions created by `npend_deriver_DAE.m` to solve a simulation.

- The function `nlink_DAE.m` (within `Lpend/`) is called by `npend.m` and calls the functions created by `nlink_deriver_DAE.m` to solve a simulation.

- The function `nmode.m` is called by `npend.m` to create linearized solutions for an n-link pendulum using the Lagrange method.

- The function `nmode_eig.m` (within `nmode/`) is called by `npend.m` and calls the functions created by `nmode.m` to find the initial conditions for normal mode oscillation.

- The script `comp.m` prompts the user to enter 2 or 3 existing simulations, which will then be animated together and compared.

- The script `comp_results.m` was written specifically to compare the times to divergence of a 10-link DAE pendulum in a specific setup.

- The script `teval.m` analyzes the performance time of a given method. The user is prompted for a solution method, number of runs to average, and list of link lengths to run.

- The script `teval_plot.m` plots the perfomance times of the solution methods.

## Descriptions

    %-------------------------------------------------------------------------%
    % npend - Jeremy Turner
    % 
    % This is the file that drives all the magic. By running npend.m, the user
    % is promped on the command line to:
    % - select to play an animation
    % --- list all available animations
    % --- play animations
    % --- replay animations
    % --- quit
    % - create a simulations
    % --- enter number of links
    % --- choose solution method
    % --- specify parameters and initial conditions
    % ----- value can be entered, or option selected to assign values
    % ------- randomized given a mean and standard deviation
    % ------- specified individually
    % --- replay animation
    % --- save animation
    % --- rerun the simulation with new parameters and initial conditions
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % animate - Jeremy Turner
    % 
    % Animates a saved simulation.
    %
    % Input: fn - filename, e.g. 'animations/DAEp100.mat'
    %
    % Loads in all the fields of the saved simulation. Animation is the same as
    % in npend.m, see there for details
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_deriver_NE - Jeremy Turner
    % 
    % Symbolic solution for the equations of motion of the n-link pendulum
    % using the Newton-Euler method (Linear and Angular momentum balances).
    %
    % Input: n - number of links
    %
    % Creates two files, npend_alphas_NE_M_[n].m and
    % npend_alphas_Lagrange_NE_[n], which are the [M] and [b] matrices such
    % that the solutions [x] = [M]^-1 [b].
    %
    % Extra: Includes simple friction force.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_deriver_L - Jeremy Turner
    % 
    % Symbolic solution for the equations of motion of the n-link pendulum
    % using the Lagrange method.
    %
    % Input: n - number of links
    %
    % Creates two files, npend_alphas_Lagrange_M_[n].m and
    % npend_alphas_Lagrange_b_[n], which are the [M] and [b] matrices such that
    % the solutions [x] = [M]^-1 [b].
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_deriver_DAE - Jeremy Turner
    % 
    % Symbolic solution for the equations of motion of the n-link pendulum
    % using the Differential Alebraic Equations method.
    %
    % Input: n = number of links
    %
    % Creates two files, npend_alphas_DAE_M_[n].m and
    % npend_alphas_Lagrange_DAE_[n], which are the [M] and [b] matrices such
    % that the solutions [x] = [M]^-1 [b].
    %
    % Extra: The base can accelerate sinusoidally in any direction, with any
    % freuency, and with any amplitude.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % nlink_deriver_DAE - Jeremy Turner
    % 
    % Symbolic solution for the equations of motion of the n-bar linkage using 
    % the Differential Alebraic Equations method.
    %
    % Within the solver, n = n - 1 such that the system is treated like a
    % pendulum with one end fixed.
    %
    % Input: n = number of links
    %
    % Creates two files, nlink_alphas_DAE_M_[n].m and
    % nlink_alphas_Lagrange_DAE_[n], which are the [M] and [b] matrices such
    % that the solutions [x] = [M]^-1 [b].
    %
    % Extra: Restoring force can be applied to the first link. It provides a
    % very small driving force in the direction of rotation (opposite of
    % friction). The intended purpose is to allow certain systems (e.g. the
    % four bar crank-rocker) to make complete rotations. This should not be
    % applied over long integration times as it continuously accelerates the
    % system.
    %
    % Extra: The end can accelerate sinusoidally in any direction, with any
    % freuency, and with any amplitude. This breaks the constraint of the nth
    % bar, so it is really a way to drive the end of an n-link pendulum.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_NE - Jeremy Turner
    % 
    % Right-hand side function called by ode45 in npend.m to simulate system
    % using the symbolically-derived ODEs using the Newton-Euler method.
    %
    % Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
    %        p - Parameter struct
    %       Mf - [M] matrix symbolically derived function file name
    %       bf - [b] vector symbolically derived function file name
    %
    % Returns: State vector at new timestep
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_Lagrange - Jeremy Turner
    % 
    % Right-hand side function called by ode45 in npend.m to simulate system
    % using the symbolically-derived ODEs using the Lagrange method.
    %
    % Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
    %        p - Parameter struct
    %       Mf - [M] matrix symbolically derived function file name
    %       bf - [b] vector symbolically derived function file name
    %
    % Returns: State vector at new timestep
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % npend_DAE - Jeremy Turner
    % 
    % Right-hand side function called by ode45 in npend.m to simulate system
    % using the symbolically-derived ODEs using the DAE method. Includes
    % shaking base.
    %
    % Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
    %        p - Parameter struct
    %        t - timespan as defined in npend.m
    %       Mf - [M] matrix symbolically derived function file name
    %       bf - [b] vector symbolically derived function file name
    %
    % Returns: State vector at new timestep, including the new position of the
    % base.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % nlink_DAE - Jeremy Turner
    % 
    % Right-hand side function called by ode45 in npend.m to simulate system
    % using the symbolically-derived ODEs using the DAE method. Includes
    % shaking end.
    %
    % Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
    %        p - Parameter struct
    %        t - timespan as defined in npend.m
    %       Mf - [M] matrix symbolically derived function file name
    %       bf - [b] vector symbolically derived function file name
    %
    % Returns: State vector at new timestep
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % nmode - Jeremy Turner
    % 
    % Linearization of symbolic solution for the equations of motion of the 
    % n-link pendulum using the Lagrange method.
    %
    % Input: n - number of links
    %
    % Creates two files, nmodeM_[n].m and nmodeK_[n].m, which are the [M] and
    % [K] matrices such that eig(M\K) returns the normal modes of oscillations
    % of the system.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % nmode_eig - Jeremy Turner
    % 
    % Function called in npend.m to find the normal modes of oscillation by
    % solving M\K using system parameters.
    %
    % Input: p - Parameter struct
    %       Mf - [M] matrix symbolically derived function file name
    %       Kf - [K] matrix symbolically derived function file name
    %
    % Returns: n x 2n matrix, where first n vectors are the initial theta
    % values for the n normal modes, and the last n columns form a diagonal
    % matrix whose elements are the squared frequencies of normal mode
    % oscillations.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % comp - Jeremy Turner
    % 
    % Compares either 2 or 3 existing simulations. The simulations must run for
    % the same amount of time / have the same number of timesteps. If 3
    % simulations are entered, they will be animated together, and their total
    % energies will be plotted together on the same figure. If 2 simulations
    % are entered, an additional plot will show when the two simulations begin
    % to diverge. This is intended for similar simulations. Divergence is
    % defined as when any pair the end points of the links between the pendula
    % become separated by more than 0.001 of the average link length.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % comp_results - Jeremy Turner
    %
    % This file was written solely to compare the time to divergence of a
    % certain set of files.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % teval - Jeremy Turner
    % 
    % This function runs the symbolic deriver methods to analyze their
    % performance. Given a number of runs to average over (m) and the numbers
    % of links to use (n), teval.m will run the specified method for n links m
    % times each. The data will be saved to a file only specifying the method
    % used, so that if teval.m is run again with the same method, the data will
    % be over-written unless the original file has been renamed. Additionally,
    % teval will plot the performance data.
    % ------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % teval_plot - Jeremy Turner
    % 
    % teval_plot.m will plot the performance data of the four methods, assuming
    % that each method has been analyzed by teval.m
    % ------------------------------------------------------------------------%
