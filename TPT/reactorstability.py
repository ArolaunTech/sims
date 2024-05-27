import random

initBREL = 0.15
reactionProg = 0.8
print(0, initBREL, 1-initBREL)

reactorTHRM = (1-2*abs(initBREL-0.5))*reactionProg
reactorBREL = initBREL - 0.5*reactorTHRM
reactorBRMT = 1-reactorTHRM-reactorBREL
print(reactorTHRM, reactorBREL, reactorBRMT)

outputBREL = reactorBREL
outputBRMT = 1-outputBREL
print(0, outputBREL, outputBRMT)

F = (1-initBREL)/outputBRMT
print(F)