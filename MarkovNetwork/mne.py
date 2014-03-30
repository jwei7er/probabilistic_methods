"""
Markov Network Evaluator

This program reads in a discrete Markove network file and prints out the partition function for the network

Author: Jordan Weiler
Date:   May 3, 2013

ex: python mne.py file.uai
"""

from factor import Factor
import sys, math

global inFile, outFile, debug

def closeFiles():
    """
    Close the input and output files
    """
    inFile.close()
    outFile.close()

def factorProduct(f1, f2, mnCards):
    """
    Multiply two factors by finding same variable assignments and multiplying their phi together.
    """
    
    # Get the unique variables across both factors 
    uniqueVars = f2.variables[:]
    for i in range(len(f1.variables)):
        if f1.variables[i] not in uniqueVars:
            uniqueVars.insert(0, f1.variables[i])

    # Create the new product factor 
    factor = Factor(uniqueVars)
    for i in range(len(factor.variables)):
        factor.setCard(i, mnCards[factor.variables[i]])
    factor.calculateStrides()

    j, k = 0, 0
    psi = []
    assignment = [0 for i in range(len(factor.variables))]
    factorSize = 1

    # Size of new factor will be product of each variable's cardinality 
    for i in range(len(factor.card)):
        factorSize *= factor.card[i]

    for i in range(factorSize):
        if debug: print "multiplying", j, k, f1.phi[j], f2.phi[k]
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
    Parses the input arguments
    """
    global inFile, outFile, debug
    debug = False

    args = sys.argv[1:]

    if len(args) < 1:
        raise Exception("Error: Input file required.")

    # Open input file
    inFile = open(args[0], "r")

    if len(args) > 1:
        for i in range(1, len(args)):
            if args[i].lower() == "debug":
                debug = True
            else:
                raise Exception("Error: " + args[i] + " argument not recognized")

    # Open output file
    outFile = open(args[0] + ".pr", "w")

def readFunctionTables(mnCards, mnCliques):
    """
    Reads function tables defined in the input file
    """
    mnFactors = []
    factorIndex, factorNum = 0, 0

    fClass = Factor(mnCliques[factorIndex])
    for i in range(len(fClass.variables)):
        fClass.setCard(i, mnCards[fClass.variables[i]])
    fClass.calculateStrides()
    
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
                        fClass = Factor(mnCliques[factorIndex])
                        for i in range(len(fClass.variables)):
                            fClass.setCard(i, mnCards[fClass.variables[i]])
                        fClass.calculateStrides()
            else:
                factorNum = int(inputRow)
                
        if len(mnFactors) >= len(mnCliques):
            # We have read in all the cliques so break out of while loop 
            done = True
    
    return mnFactors

def readInputFile():
    """
    Reads the preamble and function tables from the input file
    """
    mnCards, mnCliques = readPreamble()
	
    mnFactors = readFunctionTables(mnCards, mnCliques)
    
    return mnCards, mnFactors

def readPreamble():
    """
    Read preamble defined in the input file
    """
    mnCards, mnCliques = None, None
    
    inputRow = inFile.readline().strip()
    if inputRow.lower() == "markov":
        # Number of variables 
        mnVars = int(inFile.readline().strip())
        if debug: print "mnVars:", mnVars
        
        # Variable cardinality 
        mnCards = []
        cardinalityAr = inFile.readline().strip().split()
        for c in cardinalityAr:
            mnCards.append(int(c))
        if debug: print "mnCards:", mnCards
    
        # Variable cliques 
        numCliques = int(inFile.readline().strip())
        mnCliques = []
        for i in range(numCliques):
            cliqueAr = inFile.readline().strip().split()
            
            cliqueList = []
            for j in range(1, len(cliqueAr)):
                cliqueList.append(int(cliqueAr[j]))
            mnCliques.append(cliqueList)
        if debug: print "mnCliques:", mnCliques, "\n"

        return mnCards, mnCliques

    else:
        raise Exception("Error: " + args[i] + " argument not recognized")

def solvePR(mnCards, mnFactors):
    """
    Solve the partition function
    """
    if debug: print "solving PR"
    outFile.write("PR\n")

    if debug:
        mnFactors[0].printF()
        mnFactors[1].printF()

    if len(mnFactors) > 1:
        # Get the first product of factors
        newF = factorProduct(mnFactors[0], mnFactors[1], mnCards)
        if debug: newF.printF()

        for i in range(1, len(mnFactors)-1):
            # Combine all other factors into total product
            newF = factorProduct(newF, mnFactors[i+1], mnCards)
            
            if debug:
                print "newF"
                newF.printF()
    else:
        newF = mnFactors[0]

    # Sum all values in the factor 
    total = 0
    for i in range(len(newF.phi)):
        total += newF.phi[i]
    Z = str(total)

    print Z
    outFile.write(Z)
    	
if __name__ == "__main__":
    """
    The main function called when mne.py is run from the command line
    """
    parseInputArguments()

    mnCards, mnFactors = readInputFile()

    solvePR(mnCards, mnFactors)

    closeFiles()

