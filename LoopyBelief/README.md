Instructions for running loopyBP.py
===================================

This is a command line program that reads in a Markov network file and prints out the univariate marginals using belief propagation. The marginal distribution over each variable is printed on a single line, following the same variable order as in the file.

To run loopyBP.py from the command line, type:

        python loopyBP.py file.uai

Optional arguments are also available. To run in debug mode, type:

        python loopyBP.py file.uai debug

To print out the execution time, type:

        python loopyBP.py file.uai time
