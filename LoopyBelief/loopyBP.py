"""
Markov Network Evaluator Using Loopy Belief Propagation

This command line program reads in a Markov network file and prints out the univariate marginals using belief propagation. The marginal distribution over each variable is printed on a single line, following the same variable order as in the file.

Author: Jordan Weiler
Date:   May 31, 2013

ex. python loopyBP.py file.uai
"""

from factor import Factor
from message import Message
import sys, math, datetime

global inFile
global mnVars, mnCards, mnCliques, mnFactors, mnMarginals
global showTime, debug

def closeFiles():
    """
    Close any opened files
    """
    inFile.close()

def createMessages():
    """
    Create the variable to factor and factor to variable messages
    """
    global vToF, fToV, dictVtoF, dictFtoV
    vToF, fToV = [], []
    for i in range(len(mnFactors)):
        vrbls = mnFactors[i].variables
        for j in range(len(vrbls)):
            vToF.append(Message(vrbls[j], i, True, mnFactors[i].card[j]))
            fToV.append(Message(vrbls[j], i, False, mnFactors[i].card[j]))
            
    dictVtoF, dictFtoV = dict(), dict()
    
    for i in range(len(vToF)):
        if debug:
            print "V to F Message"
            vToF[i].printM()

        key = str(vToF[i].fact)
        if key not in dictVtoF:
            dictVtoF[key] = []
        dictVtoF[key].append(i)
    
        if debug:
            print "F to V Message"
            fToV[i].printM()

        key = str(fToV[i].var)
        if key not in dictFtoV:
            dictFtoV[key] = []
        dictFtoV[key].append(i)
    
    if debug: 
        print ""
        print "V to F"
        print dictVtoF
        print ""
        print "F to V"
        print dictFtoV
        print ""
        
def multiplyFactorAndMessages(factor, msgList, varToIgnore):
    """
    Multiply factor and messages together
    """
    if debug:
        print "Before multiplying messages"
        factor.printF()
        print "Message List"
        print msgList
        print "Ignoring ", varToIgnore
    
    for msg in msgList:
        if vToF[msg].var != varToIgnore:
            varIndex = factor.variables.index(vToF[msg].var)
            
            updatedPhi = [False for x in range(factor.size)]
            
            numUpdates = factor.size / factor.card[varIndex]
            
            for i in range(factor.card[varIndex]):
                phiIndex = 0
                for j in range(factor.size):
                    if updatedPhi[j] == False:
                        phiIndex = j
                        break
                
                updates = 0
                while updates < numUpdates:
                    strideCnt = 0
                    while strideCnt < factor.stride[varIndex]:
                        factor.phi[phiIndex] *= vToF[msg].val[i]
                        updatedPhi[phiIndex] = True
                        
                        strideCnt += 1
                        updates += 1
                        phiIndex += 1
                    phiIndex += factor.stride[varIndex] * (factor.card[varIndex] - 1)
    if debug: 
        print "After multiplying messages"
        factor.printF()
    
    return factor
        
def parseInputArguments():
    """
    Read in all input arguments and set global variables
    """
    global inFile, showTime, debug
    showTime = False
    debug = False

    args = sys.argv[1:]
    
    if len(args) < 1:
        raise Exception("Error: Input file required.")

    inFile = open(args[0], "r")

    if len(args) > 1:
        for i in range(1, len(args)):
            if args[i].lower() == "time":
                showTime = True
            elif args[i].lower() == "debug":
                debug = True
            else:
                raise Exception("Error: " + args[i] + " argument not recognized")

def printMarginals():
    """
    Print out the marginals of the Markov network
    """
    marginals = ""
    
    for i in range(mnVars):
        msg = dictFtoV[str(i)]
        
        varProb = []
        varTot = 0.0
        size = len(msg)
        if size > 0:
            card = len(fToV[msg[0]].val)
            
            for j in range(card):
                tot = 1.0
                for k in range(size):
                    tot *= fToV[msg[k]].val[j]
                varProb.append(tot)
                varTot += tot

        for j in range(len(varProb)):
            varProb[j] = varProb[j] / varTot
        
        for x in varProb:
            marginals += str(x) + " "
        marginals += "\n"
    
    print marginals 

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
                    # Reset for next factor
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
            # All cliques have been read so break out of while loop 
            done = True
    if debug:
        for i in range(len(mnFactors)):
            mnFactors[i].printF()

def readInputFile():
    """
    Read the preamble and function tables from the input file
    """
    readPreamble()
	
    readFunctionTables()

def readPreamble():
    """
    Reads preamble defined in the input file
    """
    global mnVars, mnCards, mnCliques
    
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
    else:
        raise Exception("Error: Markov network input file format needed")

def solveLoopyBP():
    """
    Solve loopy belief propagation
    """
    createMessages()
    
    loop = 0
    done = False
    while loop < 50 and not done:
        # Save current F to V and V to F lists
        sFtoV, sVtoF = [], []
        for i in range(len(fToV)):
            for j in range(len(fToV[i].val)):
                sFtoV.append(fToV[i].val[j])
            for j in range(len(vToF[i].val)):
                sVtoF.append(vToF[i].val[j])
        
        if debug: 
            print sFtoV
            print sVtoF
        
        # Update messages in both directions
        updateFtoVMessages()
        updateVtoFMessages()
       
        # Compare differences with last marginals
        diff = 0.0
        index = 0
        for i in range(len(fToV)):
            for j in range(len(fToV[i].val)):
                diff += math.fabs(fToV[i].val[j] - sFtoV[index])
                index +=1
        index = 0
        for i in range(len(vToF)):
            for j in range(len(vToF[i].val)):
                diff += math.fabs(vToF[i].val[j] - sVtoF[index])
                index += 1
        
        if debug: print "diff:", diff
            
        # If differences fall below 0.00001 stop
        if diff < 0.00001:
            done = True
            
        loop += 1
    
    if debug: print "Loops: ", loop, "\n"
    
    printMarginals()

def sumOutVariable(factor, variable):
    """
    Sum out a variable from a factor
    """
    newVars = [x for x in factor.variables if x != variable]
    varIndex = factor.variables.index(variable)
    newF = Factor(newVars, mnCards)
    
    if debug: newF.printF()
    
    usedVar = [False for x in range(factor.size)]
    
    for i in range(newF.size):
        newF.phi.append(0)
        
        start = 0
        for k in range(factor.size):
            if usedVar[k] == False:
                start = k
                break
        
        for j in range(factor.card[varIndex]):
            key = start + factor.stride[varIndex] * j
            newF.phi[i] += factor.phi[key]
            usedVar[key] = True

        if debug: print "phi: ", i, " ", newF.phi[i]

    return newF
    
def updateFtoVMessages():
    """
    Update factor to variable messages
    """
    global fToV
    
    if debug: print "## Updating Factor to Variable Message"
    
    for i in range(len(fToV)):
        if debug: 
            fToV[i].printM()

        factor = mnFactors[fToV[i].fact]
        
        msgs = dictVtoF[str(fToV[i].fact)]
        if debug:
            print "msgs: ", msgs
        
        factor = multiplyFactorAndMessages(factor, msgs, fToV[i].var)

        factorVarLen = len(factor.variables)        
        if factorVarLen > 1:
            while factorVarLen > 1:
                for j in range(factorVarLen):
                    if factor.variables[j] != fToV[i].var:
                        factor = sumOutVariable(factor, factor.variables[j])
                        factorVarLen = len(factor.variables)
                        break
        factor.renormalize()
        if debug:
            print ""
            print "After Factor Normalization"
            factor.printF()
            print "var: ", fToV[i].var, "fact: ", fToV[i].fact
        
        fToV[i].val = factor.phi[:]
        fToV[i].renormalize()
        
        if debug:
            print "#Message after normalization"
            fToV[i].printM()
    if debug: print ""

def updateVtoFMessages():
    """
    Update variable to factor messages
    """
    global vToF
    for i in range(len(vToF)):
        msgs = dictFtoV[str(vToF[i].var)]
        msgToMult = []
        for j in range(len(msgs)):
            if vToF[msgs[j]].fact != vToF[i].fact:
                msgToMult.append(vToF[msgs[j]])
        
        if len(msgToMult) > 0:
            size = len(msgToMult)
            card = len(msgToMult[0].val)
            
            for j in range(card):
                tot = 1.0
            
                for k in range(size):
                    tot *= msgToMult[k].val[j]
                vToF[i].val[j] = tot
        
    for i in range(len(vToF)):
        vToF[i].renormalize()
    
if __name__ == "__main__":
    """
    The main function called when ve.py is run from the command line
    """
    parseInputArguments()

    readInputFile()

    startTime = datetime.datetime.now()
    
    solveLoopyBP()

    endTime = datetime.datetime.now()
    
    if showTime: print "execution time:", (endTime - startTime)
       
    closeFiles()

