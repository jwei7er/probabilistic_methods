"""
Factor class holds variables, cardinality, stride values, and phi table for a factor.

Author: Jordan Weiler
Date:   May 31, 2013
"""

class Factor:
    def __init__(self, variables, cardValues):
        self.variables = variables
        self.card = []
        self.stride = []
        self.phi = []
        self.size = 1
        self.setCards(cardValues)
        self.calculateStrides()

    def calculateStrides(self):
        """
        Calculate the strides for each of the variables and the total size of the factor
        """
        for i in range(len(self.card)-1, -1, -1):
            if i == len(self.card) - 1:
                self.stride[i] = 1
            else:
                total = 1
                for j in range(i+1, len(self.card)):
                    total *= self.card[j]
                self.stride[i] = total
            self.size *= self.card[i]

    def printF(self):
        """
        Print the factor details
        """
        print "Factor"
        print "var:", self.variables
        print "card:", self.card
        print "stride:", self.stride
        print "size:", self.size
        print "phi:", self.phi
        print ""

    def setCard(self, cardValue):
        """
        Set the cardinality for a variable
        """
        self.card.append(cardValue)
        self.stride.append(0)

    def setCards(self, cardValues):
        """
        Sets the cardinality for each of the variables
        """
        for i in range(len(self.variables)):
            self.setCard(cardValues[self.variables[i]])
    
    def renormalize(self):
        """
        Renormalize the phi value
        """
        tot = 0.0
        for i in range(len(self.phi)):
            tot += self.phi[i]
        for i in range(len(self.phi)):
            self.phi[i] = self.phi[i] / tot

