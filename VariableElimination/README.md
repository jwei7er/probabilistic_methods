Instructions for running ve.py
==============================

This is a command line program that reads in a discrete Markov network file and prints out the partition function for the network. The partition function is computed by multiplying all factors together and adding up all the resulting values.

1. To run mne.py from the command line, type:

        python ve.py file.uai

2. Optional arguments can be provided after the input file (order and case do not matter).

    a. *time* - This will show you the execution time.

    b. *order=#* - Replace # with a valid ordering heuristic. This allows you to select the ordering heuristic used in variable elimination. Min-neighbors is used by default.
    Valid ordering heuristics are:
        Neighbor
        Weight
        Fill
        WgtFill

    Example of use: *order=Neighbor*

    c. *debug* - This will show you print statements useful when debugging.

    An example command line using all optional arguments is:

        python ve.py file.uai time order=fill debug

