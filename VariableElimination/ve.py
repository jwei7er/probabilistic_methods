"""
Markov Network Evaluator Using Variable Elimination

This command line program reads in a Markov network file and prints out the partition function for the network (computed with variable elimination). The variable elimination algorithm uses the min-neighbors heuristic for variable ordering by default.

Author: Jordan Weiler
Date:   May 17, 2013

ex. python ve.py file.uai
"""

from factor import Factor
import sys, math

global inFile
global mnVars, mnCards, mnCliques, mnFactors
global debug

def closeFiles():
    """
    Close the input file
    """
    inFile.close()

def factorProduct(f1, f2):
    """
    Multiplies two factors by finding same variable assignments and multiplying their phi together
    """
    
    # Get the unique variables across both factors
    uniqueVars = f2.variables[:]
    for i in range(len(f1.variables)):
        if f1.variables[i] not in uniqueVars:
            uniqueVars.insert(0, f1.variables[i])

    # Create the new product factor
    factor = Factor(uniqueVars, mnCards)

    j, k = 0, 0
    psi = []
    assignment = [0 for i in range(len(factor.variables))]
    
    for i in range(factor.size):
        #if debug: print "multiplying", j, k, f1.phi[j], f2.phi[k]
        psi.append(f1.phi[j] * f2.phi[k])
        
        # Loop from last to first since values are 000, 001, 010, 011 instead of 000, 100, 010, 110
        for v in range(len(factor.variables)-1,-1,-1):
            curVariable = factor.variables[v]
            assignment[v] += 1
            if assignment[v] == factor.card[v]:
                assignment[v] = 0
                if curVariable in f1.variables:
                    j = j - (factor.card[v] - 1) * f1.stride[f1.variables.index(curVariable)]
                if curVariable in f2.variables:
                    k = k - (factor.card[v] - 1) * f2.stride[f2.variables.index(curVariable)]
            else:
                if curVariable in f1.variables:
                    j = j + f1.stride[f1.variables.index(curVariable)]
                if curVariable in f2.variables:
                    k = k + f2.stride[f2.variables.index(curVariable)]
                break
                
    factor.phi = psi[:]

    return factor
    
def parseInputArguments():
    """
    Parse the input arguments
    """
    global inFile, debug
    debug = False

    args = sys.argv[1:]

    if len(args) < 1:
        raise Exception("Error: Input file required.")

    inFile = open(args[0], "r")

    if len(args) > 1:
        for i in range(1, len(args)):
            if args[i].lower() == "debug":
                debug = True
            else:
                raise Exception("Error: " + args[i] + " argument not recognized")

def readFunctionTables():
    """
    Reads function tables defined in the input file
    """
    global mnFactors
    mnFactors = []
    factorIndex, factorNum = 0, 0

    fClass = Factor(mnCliques[factorIndex], mnCards)
    
    factor = []
    done = False
    while not done:
        inputRow = inFile.readline().strip()
        if inputRow == "":
            continue
        else:
            if " " in inputRow:
                splitRow = inputRow.split()
                if factorNum == 0:
                    factorNum = int(splitRow[0])
                    for i in range(1, len(splitRow)):
                        factor.append(float(splitRow[i]))
                else:
                    for i in range(len(splitRow)):
                        factor.append(float(splitRow[i]))
                        
                if factorNum == len(factor):
                    factorNum = 0
                    fClass.phi = factor
                    mnFactors.append(fClass)
                    factor = []
                    factorIndex += 1
                    if factorIndex < len(mnCliques):
                        fClass = Factor(mnCliques[factorIndex], mnCards)
            else:
                factorNum = int(inputRow)
                
        if len(mnFactors) >= len(mnCliques):
            """ read in all the cliques so break out of while loop """
            done = True

def readInputFile():
    """
    Read the preamble and function tables from the input file
    """
    readPreamble()
    	
    readFunctionTables()

def readPreamble():
    """
    Read preamble defined in the input file
    """
    global mnVars, mnCards, mnCliques
    
    inputRow = inFile.readline().strip()
    if inputRow.lower() == "markov":
        """ number of variables """
        mnVars = int(inFile.readline().strip())

        if debug: print "mnVars:", mnVars
        
        """ variable cardinality """
        mnCards = []
        cardinalityAr = inFile.readline().strip().split()
        for c in cardinalityAr:
            mnCards.append(int(c))

        if debug: print "mnCards:", mnCards
    
        """ variable cliques """
        numCliques = int(inFile.readline().strip())
        mnCliques = []
        for i in range(numCliques):
            cliqueAr = inFile.readline().strip().split()
            
            cliqueList = []
            for j in range(1, len(cliqueAr)):
                cliqueList.append(int(cliqueAr[j]))
            mnCliques.append(cliqueList)

        if debug: print "mnCliques:", mnCliques, "\n"
    else:
        raise Exception("Error: Markov network input file needed")

def factorProductList(factors):
    """
    Find the product of a list of factors
    """
    if debug: print "factor product list: ", len(factors)

    if len(factors) < 2:
        return factors[0]
    else:
        newF = factors[0]
        for i in range(1, len(factors)):
            newF = factorProduct(newF, factors[i])
        return newF
    
def sumOutVariable(factor, variable):
    """
    Sum out a variable for variable elimination
    """
    if debug: print "eliminating factor"

    newVars = [x for x in factor.variables if x != variable]

    if debug: print "new vars: ", newVars

    newF = Factor(newVars, mnCards)

    if len(newVars) > 0:
        varIndex = factor.variables.index(variable)
        if debug: 
            print "var index: ", varIndex
            newF.printF()

        usedVar = [False for x in range(factor.size)]

        if debug: 
            print "used var: ", usedVar
            print "stride: ", factor.stride[varIndex], ", card: ", factor.card[varIndex]

        psi = []
        for i in range(newF.size):
            psi.append(0)
    
            start = 0
            for k in range(factor.size):
                if usedVar[k] == False:
                    start = k
                    break
    
            for j in range(factor.card[varIndex]):
                if debug: print "start: ", start, " stride: ", factor.stride[varIndex], " j: ", j
                psi[i] += factor.phi[start + factor.stride[varIndex] * j]
                usedVar[start + factor.stride[varIndex] * j] = True

            if debug: print "psi: ", i, " ", psi[i]

        newF.phi = psi[:]

    return newF
    
def solvePR():
    """
    Solve the partition function
    """
    if debug: print "solving PR"

    if len(mnFactors) > 1:
        for i in range(mnVars - 1):
            if debug: print "eliminating : ", i

            eliminateV = i
            elimIndexSet = []
            elimFactors = []

            for j in range(len(mnFactors)-1, -1, -1):
                if debug: print "j: ", mnFactors[j].variables
                if eliminateV in mnFactors[j].variables:
                    if debug: print "adding : ", j
                    elimIndexSet.append(j)
                    elimFactors.append(mnFactors[j])
                    mnFactors.pop(j)
                    
            if debug: print "elim Set: ", elimIndexSet

            newF = factorProductList(elimFactors)

            if debug:
                print "newF: "
                newF.printF()
            
            newF = sumOutVariable(newF, eliminateV)

            if debug:
                print "improvedF: "
                newF.printF()

            mnFactors.append(newF)

        newF = mnFactors[0]

        if debug: 
            print "Last factor:"
            newF.printF()
    else:
        newF = mnFactors[0]

    # Sum all values in the factor
    total = 0
    for i in range(len(newF.phi)):
        total += newF.phi[i]
    Z = str(total)

    print Z
    	
if __name__ == "__main__":
    """
    The main function called when ve.py is run from the command line
    """
    parseInputArguments()

    readInputFile()

    solvePR()

    closeFiles()

