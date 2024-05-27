# sims
My simulations and optimizers. Please note that these were built for myself, not for anyone else, so it might be hard to use. No tech support will be provided by me.

A lot of these programs are unfinished and never will be completed.

## Brief Code Descriptions

### KSP simulations:
 - ```comms.py```: Optimizes satellite communication networks with a given number of satellites for maximum reliability.
 - ```fuelopt.py```: Optimizes the distribution of fuel across rocket stages for maximum delta-v.
 - ```graphtraverse.py```: Optimizes the location of a mothership in the Jool system for the most efficient Jool-5. (yes I know it's Elliptical Laythe and Tylo Orbit)
 - ```ion.go```: Optimizes a Moho SSTO for minimum mass.
 - ```ion_engine.py```: Optimizes ion chairs for minimum mass while still being able to go to Mun.
 - ```j5ssto.py```: Optimizes Jool-5 SSTO route for maximum fractional payload.
 - ```planesim.cpp```: Optimizes RAPIER-based ascent profiles to reach a given suborbital trajectory using minimal fuel. One of my first C++ programs so may contain bad practices.
 - ```tylo_landing.py```: Optimizes Tylo landers for minimum mass.

### Minecraft simulation:
 - ```Wheat.py```: Finds optimal time to harvest a wheat farm. One of my first programs so isn't very good.

### TPT simulations:
 - ```7seg.py```: Finds wiring designs for 7-segment displays.
 - ```countercurrent.py```: Simulates countercurrent heat exchangers.
 - ```draybomb.py```: Simulates a DRAY bomb.
 - ```fuse.py```: Optimizes platinum-catalyzed cold fusion reaction rate.
 - ```logic.py```: Optimizes arrangements of RSST-based logic gates that perform certain functions for least complexity. Used to design my RSST-based XOR gate and half adder.
 - ```reactormodel.py```, ```reactorstability.py```: Simulates thermite reactors and shows the space (thankfully 2-dimensional!) of all possible THRM reactors.
 - ```vibr.py```: Simulates VIBR spark generation.