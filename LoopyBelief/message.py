"""
Message class holds a message between a factor and a variable

Author: Jordan Weiler
Date:   May 31, 2013
"""

class Message:
    def __init__(self, var, fact, toFact, size):
        self.var = var          #variable
        self.fact = fact        #factor
        self.toFact = toFact
        if self.toFact:
            self.f = self.var   #from
            self.t = self.fact  #to
        else:
            self.f = self.fact  #from
            self.t = self.var   #to
        self.val = []
        for i in range(size):
            self.val.append(1.0)

    def renormalize(self):
        """
        Renormalize the values to add up to 1.0
        """
        tot = 0.0
        for i in range(len(self.val)):
            tot += self.val[i]

        for i in range(len(self.val)):
            if tot == 0:
                self.val[i] = 1.0 / len(self.val)
            else:
                self.val[i] = self.val[i] / tot

    def printM(self):
        """
        Print out the message details
        """
        if self.toFact:
            print "M " + str(self.f) + " -> F" + str(self.t) + " (" + str(self.var) + ") = " + str(self.val)
        else:
            print "M F" + str(self.f) + " -> " + str(self.t) + " (" + str(self.var) + ") = " + str(self.val)

