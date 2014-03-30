"""
Factor class holds variables, cardinality, stride values, and phi table for a factor.

Author: Jordan Weiler
Date:   May 3, 2013
"""

class Factor:
    def __init__(self, variables):
        self.variables = variables
        self.card = []
        self.stride = []
        self.phi = []

    def calculateStrides(self):
        """
        Calculate the strides for each of the variables
        """
        for i in range(len(self.card)-1, -1, -1):
            if i == len(self.card) - 1:
                self.stride[i] = 1
            else:
                total = 1
                for j in range(i+1, len(self.card)):
                    total *= self.card[j]
                self.stride[i] = total

    def printF(self):
        """
        Print the factor details
        """
        print "Factor"
        print "var:", self.variables
        print "card:", self.card
        print "stride:", self.stride
        print "phi:", self.phi
        print ""

    def setCard(self, cardIndex, cardValue):
        """
        Set the cardinality at a specific variable index
        """
        self.card.append(cardValue)
        self.stride.append(0)


