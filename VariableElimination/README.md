Instructions for running ve.py
==============================

This is a command line program that reads in a Markov network file and prints out the partition function for the network (computed with variable elimination). The variable elimination algorithm uses the min-neighbors heuristic for variable ordering.

To run ve.py from the command line, type:

        python ve.py file.uai

An optional argument for debugging is also available. To run in debug mode, type:

        python ve.py file.uai debug

