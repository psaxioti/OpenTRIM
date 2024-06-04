- Fix the seed for the rng \
  1 seed value (an int) is needed to set the state of the rng. In case of multiple threads, use rng.long_jump to set the state in each thread

- Save the rng state to the h5 output file

- Continue a previous simulation by reading the h5 output file. Stream files should append new events

- Make user-defined tallies for various events. The tallies should have the possibility to have their own mesh 

- Handle surface effects