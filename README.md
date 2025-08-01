# sims
Some simulations and optimizers. Please note that these were built for myself, not for anyone else, so it might be hard to use. No tech support will be provided by me.

A lot of these programs are unfinished and never will be completed.

## Brief Code Descriptions

### Real life simulations:
 - ```1d.py```: An attempt at simulating the inside of a star.
 - ```bloom.py```: Simulates the appearance of bright objects with different cameras/eyes.
 - ```brownian.py```: Simulates hydraulic erosion and (badly) simulates climate. Not realistic at all but it works.
 - ```cmt.py```: Shows different colormaps so they can be compared.
 - ```collide.py```: A highly abstracted simulation of molecular collisions.
 - ```convect.py```: An attempt at simulating convection.
 - ```detail.py```: An unfinished raytracer.
 - ```diffuse.py```: An attempt at simulating fluids.
 - ```equilibrium.py```: See ```collide.py```.
 - ```fluid.py```: Yet another attempt at simulating fluids.
 - ```graphs.py```: Finds day lengths at different latitudes and times of year.
 - ```hum.py```: Imitates the sound of some weird string instrument.
 - ```lbm.py```, ```lbm2.py```: Two successful attempts at simulating fluids using the Lattice Boltzmann method.
 - ```n-body.py```: Simulates n-body gravitation.
 - ```panel.py```: Simulates aerodynamics using the panel method.
 - ```plantae.py```: Simulates the inner workings of plants to determine their colors. Works okay but not the best.
 - ```render.py```: See ```detail.py```.

### KSP simulations:
 - ```autogyro.py```, ```autogyro2.py```: Calculates the aerodynamic properties of autogyros.
 - ```comms.py```: Optimizes satellite communication networks with a given number of satellites for maximum reliability.
 - ```consts.h```: Some constants. Required by some of the C++ optimizers.
 - ```curves.h```: Defines some curves. Required by some of the C++ optimizers.
 - ```fuelopt.py```: Optimizes the distribution of fuel across rocket stages for maximum delta-v.
 - ```generallandingprofile.cpp```: (Tries to) find the optimal landing profile for a given craft to land on a given body.
 - ```graphtraverse.py```: Optimizes the location of a mothership in the Jool system for the most efficient Jool-5. (yes I know it's Elliptical Laythe and Tylo Orbit)
 - ```inclination.py```: Calculates the delta-v required to change the inclination of an orbit.
 - ```ion.go```: Optimizes a Moho SSTO for minimum mass.
 - ```ion_engine.py```: Optimizes ion chairs for minimum mass while still being able to go to Mun.
 - ```ion_landing.py```: Simulates the Constant Altitude and Constant Horizontal Acceleration landing profiles.
 - ```ion_maximizer.py```: Maximizes the thrust of ion engines undergoing EC starvation.
 - ```j5ssto.py```: Optimizes Jool-5 SSTO route for maximum fractional payload.
 - ```kicks.py```: Optimizes periapsis kicks for closest encounter.
 - ```landingprofileopt.cpp```: Given fuel amounts in orbit, optimizes landing profile to use as little fuel as possible.
 - ```landingtests.py```: Simulates the fuel usage from different landing trajectories.
 - ```landingtestsgeneral.py```: Simulates the fuel usage from different landing trajectories.
 - ```planesim.cpp```: Optimizes RAPIER-based ascent profiles to reach a given suborbital trajectory using minimal fuel. One of my first C++ programs so may contain bad practices.
 - ```prop.py```: (Unsuccessfully) calculates aspects of KSP helicopters.
 - ```rapiermaxspeed.py```: (Somewhat successfully) calculates parameters of RAPIER-based Kerbin SSTOs.
 - ```relvelmagnifier.py```: Calculates gravity assist trajectories to save dV through orbit rotation shenanigans.
 - ```shift_assist.py```: Part of a multiple gravity assist planner.
 - ```singlestagedv.py```: Calculates optimal engine configuration for a single stage craft with certain TWR and vacuum delta-V. Assumes continuous fuel tanks.
 - ```single_stage.py```: A class that represents a single stage craft.
 - ```single_stage_simulation.py```: Optimizes single stage craft landings in vacuum for maximum payload fraction.
 - ```twrlanding.py```: Calculates losses in vacuum landings based on the TWR of craft.
 - ```tylo_landing.py```: Optimizes Tylo landers for minimum mass.
 - ```tylo_optimizer_test.py```: An example of how to use ```single_stage_simulation.py``` to optimize a Tylo lander.
 - ```utils.py```: A utils library for ```single_stage.py```.
 - ```vacuumoptimizer.cpp```: Optimizes single-stage vacuum landers for maximum fractional payload by changing engine count and throttle.
 - ```veltest.py```: Optimizes landing trajectories to use as little fuel as possible.

### Minecraft simulations:
 - ```Wheat.py```: Finds optimal time to harvest a wheat farm. One of my first programs so isn't very good.

### TPT simulations:
 - ```7seg.py```: Finds wiring designs for 7-segment displays.
 - ```countercurrent.py```: Simulates countercurrent heat exchangers.
 - ```draybomb.py```: Simulates a DRAY bomb.
 - ```fuse.py```: Optimizes platinum-catalyzed cold fusion reaction rate.
 - ```gravitydescent.py```: Tries to create a given force field with gravity sources.
 - ```heat.py```: Calculates the best wall setup to resist overheating.
 - ```logic.py```: Optimizes arrangements of RSST-based logic gates that perform certain functions for least complexity. Used to design my RSST-based XOR gate and half adder.
 - ```reactormodel.py```, ```reactorstability.py```: Simulates thermite reactors and shows the space (thankfully 2-dimensional!) of all possible THRM reactors.
 - ```sltw.py```: Maximizes the output of water generators.
 - ```vibr.py```: Simulates VIBR spark generation.

### Other simulations:
 - ```reaction-diffusion.py```: An attempt at simulating reaction-diffusion systems.
