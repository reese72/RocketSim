# RocketSim
The name is pretty self explanatory, this is a model rocket simulation written in python, it's not perfectly accurate but it can rival OpenRocket

# What is This?
I guess I should probably give a quick explanation as to what this simulation is, what it does, what it can do, and what it can't do. As mentioned in the repo's 
description this is a model rocket simulation that is written in python. I takes custom inputs for almost every parameter of the rocket but is not visual like 
OpenRocket. This simulation's best use case is in the form of an API of sorts, OpenRocket is objectively the more accurate simulation, but this is made in 
perfectly readable and rework-able python and can therefore easily be injected into other codebases. It's accuracy is a bit iffy in it's current state, it can 
accurately simulate rockets with up to a D-class motor, but anything above that begins to become inaccurate, I'm working on finding the bug but it is elusive. 
The simulation also makes some assumptions, it assumes a constant burn rate, drag coefficent, and cross-sectional area, which may or may not be a deal-breaker when simulation rockets.
