# Module TermSelection.py processing
# -*- coding: utf-8 -*-

import math, codecs, string
import fnmatch


print "Read TermSelection.py"


   def doAll(aFileM, aFileF, aYear="2013"):
   aFileName = "output" + aYear
   
   aL = doDifferenceSymetric("CHI", aFileM,  aFileF, aFileName+"CHI.txt", 550, 0)
   aL = doDifferenceSymetric("GSS", aFileM,  aFileF, aFileName+"GSS.txt", 550, 0)
   aL = doDifferenceSymetric("CC", aFileM, aFileF, aFileName+"CC.txt", 550, 0);
   aL = doDifferenceSymetric("GR", aFileM, aFileF, aFileName+"GR.txt", 550, 0);
   aL = doDifferenceSymetric("IG", aFileM, aFileF, aFileName+"IG.txt", 550, 0);
   
   aL = doDifferenceAsymetric("RS", aFileM, aFileF, aFileName+"RS.", 550, 0);
   aL = doDifferenceAsymetric("PMI", aFileM, aFileF, aFileName+"PMI.", 550, 0);
   aL = doDifferenceAsymetric("WLLR", aFileM, aFileF, aFileName+"WLLR.", 550, 0);
   aL = doDifferenceAsymetric("RF", aFileM,  aFileF, aFileName+"RF.", 550, 0)
   aL = doDifferenceAsymetric("OR", aFileM,  aFileF, aFileName+"OR.", 550, 0)
   return(True)


#
# Difference between two dictionaries
#

def doDifferenceProbD(aFileNameBase, aFileNameComp, outputFile="output.txt", aNumber=100, minDF=5):
   (aDictM, aNm) = createDictionary(aFileNameBase)
   print len(aDictM),
   if (minDF > 0):
      filterDictTF(aDictM, minDF)
      print len(aDictM)
   else:
      print
   aDictM = normalizeSimpleDict(aDictM)
   (aDictF, aNf) = createDictionary(aFileNameComp)
   print len(aDictF),
   if (minDF > 0):
      filterDictTF(aDictF, minDF)
      print len(aDictF)
   else:
      print
   aDictF = normalizeSimpleDict(aDictF)
   
   aD     = diffDictionary(aDictM, aDictF)
   if (aNumber > 0):
      aT1 = extractTopSortedDict(aD, aNumber, True)
      dumpDictionary(aT1, outputFile+"1M.txt")
   else:
      aT1 = aD
      dumpDictionaryAll(aD, outputFile+"1M.txt", True)
      
   aD = diffDictionary(aDictF, aDictM)
   if (aNumber > 0):
      aT2 = extractTopSortedDict(aD, aNumber, True)
      dumpDictionary(aT2, outputFile+"2F.txt")
   else:
      aT2 = aD
      dumpDictionaryAll(aD, outputFile+"2F.txt", True)      
   return(aT1, aT2)

#
# Difference with a given method
#
def doDifferenceAsymetric(aMethod, aFileNameBase, aFileNameComp, outputFile="outputOR.txt", aNumber=100, minDF=5):
   (aDictM, aNm) = createDFDictionary(aFileNameBase)
   print len(aDictM),
   if (minDF > 0):
      filterDictTF(aDictM, minDF)
      print len(aDictM)
   else:
      print
   (aDictF, aNf) = createDFDictionary(aFileNameComp)
   print len(aDictF),
   if (minDF > 0):
      filterDictTF(aDictF, minDF)
      print len(aDictF)
   else:
      print
   if (aMethod == "OR"):
      aD = computeOddRatio(aDictM, aDictF, aNm, aNf)
   if (aMethod == "PMI"):
      aD = computePMI(aDictM, aDictF, aNm, aNf)
   if (aMethod == "RF"):
      aD = computeRF(aDictM, aDictF, aNm, aNf)
   if (aMethod == "RS"):
      aD = computeRS(aDictM, aDictF, aNm, aNf)
   if (aMethod == "WLLR"):
      aD = computeWLLR(aDictM, aDictF, aNm, aNf)
      
   aT1 = extractTopSortedDict(aD, aNumber)
   dumpDictionary(aT1, outputFile+"1M.txt")
   if (aMethod == "OR"):
      aD = computeOddRatio(aDictF, aDictM, aNf, aNm)
   if (aMethod == "PMI"):
      aD = computePMI(aDictF, aDictM, aNf, aNm)
   if (aMethod == "RF"):
      aD = computeRF(aDictF, aDictM, aNf, aNm)
   if (aMethod == "RS"):
      aD = computeRS(aDictF, aDictM, aNf, aNm)
   if (aMethod == "WLLR"):
      aD = computeWLLR(aDictF, aDictM, aNf, aNm)
      
   aT2 = extractTopSortedDict(aD, aNumber)
   dumpDictionary(aT2, outputFile+"2F.txt")      
   return(aT1, aT2)

def doDifferenceSymetric(aMethod, aFileNameBase, aFileNameComp, outputFile="output.txt", aNumber=100, minDF=5):
   (aDictM, aNm) = createDFDictionary(aFileNameBase)
   print len(aDictM),
   if (minDF > 0):
      filterDictTF(aDictM, minDF)
      print len(aDictM)
   else:
      print
   (aDictF, aNf) = createDFDictionary(aFileNameComp)
   print len(aDictF),
   if (minDF > 0):
      filterDictTF(aDictF, minDF)
      print len(aDictF)
   else:
      print
   if (aMethod == "GSS"):
      aD = computeGSS(aDictM, aDictF, aNm, aNf)
   if (aMethod == "CHI"):
      aD = computeCHI(aDictM, aDictF, aNm, aNf)
   if (aMethod == "GR"):
      aD = computeGR(aDictM, aDictF, aNm, aNf)
   if (aMethod == "IG"):
      aD = computeIG(aDictM, aDictF, aNm, aNf)
   if (aMethod == "CC"):
      aD = computeCC(aDictM, aDictF, aNm, aNf)
   aDictF = aDictM = 0
   aB = extractBottomSortedDict(aD, aNumber)
   aT = extractTopSortedDict(aD, aNumber)
   dumpDictionary(aT, outputFile)
   dumpDictionary(aB, outputFile)
   return(aT, aB)

#
# Term Selection procedure 
#
#
# GSS  return a Dict (word:gss) from a corpusDFDict, subcorpusDFDict, 
def computeRF(aDict1DF, aDict2DF, subCorpusSize1, subCorpusSize2, aSet=[]):
   aDictValue = {}
   acValue = float(subCorpusSize1)
   bdValue = float(subCorpusSize2)
   nValue = acValue + bdValue
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0
      maxValue = max(bValue,1.0)
#     print aKey, aValue, bValue, maxValue, math.log (2 + (aValue/maxValue),2)
      aDictValue [aKey] =  math.log (2 + (aValue/maxValue),2)
   return(aDictValue)

def computeGSS(aDict1DF, aDict2DF, subCorpusSize1, subCorpusSize2, aSet=[]):
   aDictValue = {}
   acValue = float(subCorpusSize1)
   bdValue = float(subCorpusSize2)
   nValue = acValue + bdValue
   n2Value = nValue * nValue
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0001
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0001
      cValue = acValue - aValue
      dValue = bdValue - bValue
#     print aKey, aValue, bValue, cValue, dValue, (((aValue*dValue) - (cValue*bValue)) / n2Value)
      aDictValue [aKey] = ((aValue*dValue) - (cValue*bValue)) / n2Value
   return(aDictValue)

def computeOddRatio(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   acValue = float(corpusSize1)
   bdValue = float(corpusSize2)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0001
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0001
      cValue = acValue - aValue
      dValue = bdValue - bValue
      prob1 = aValue / acValue
      prob2 = bValue / bdValue
      if ((bValue <= 0) or (abs(aValue - acValue) <= 0.00000001)) :
         oddValue = 0.0
      else:
         oddValue = (prob1*(1-prob2)) / ((1-prob1)*prob2)
#      print aKey, aValue, bValue, cValue, dValue, oddValue
      aDictValue[aKey] = oddValue
   return(aDictValue)

def computeOddRatioOLD(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0001
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0001
      cValue = float(corpusSize1) - aValue
      dValue = float(corpusSize2) - bValue
      if ((cValue*bValue) <= 0):
         oddValue = 0.0
      else:
         oddValue = (aValue*dValue) / (cValue*bValue)
#      print aKey, aValue, bValue, cValue, dValue, oddValue
      aDictValue[aKey] = oddValue
   return(aDictValue)

# CHI Chi-square value  return a Dict (word:chi) from a corpusDFDict, subcorpusDFDict, 
def computeCHI(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue = float(corpusSize1) + corpusSize2
   acValue = float(corpusSize1)
   bdValue = float(corpusSize2)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0
      abValue = aValue + bValue
      cValue = acValue - aValue
      dValue = bdValue - bValue
      num = nValue * pow( (float(aValue) * dValue)  - (cValue * bValue), 2)
      denum = float(abValue) * acValue * (bValue+dValue) * (cValue + dValue)
      if ((num <= 0) | (denum <= 0)):
         chiValue = 0.0
      else:
         chiValue =   num / denum
      aDictValue [aKey] = chiValue
#      print aKey, aValue, bValue, cValue, dValue, chiValue
   return(aDictValue)

# PMI Point Mutual Information   return a Dict (word:pmi) from a corpusDFDict, subcorpusDFDict, 
def computePMI(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue = float(corpusSize1 + corpusSize2)
   acValue = float(corpusSize1)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0001
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0001
      abValue = aValue + bValue
      if (((abValue * acValue) <= 0) | (aValue <= 0)):
         pmiValue = 0.0
      else:
         pmiValue = math.log (  (float (aValue*nValue) / (abValue * acValue)), 2)
      aDictValue[aKey] = pmiValue
#      print aKey, aValue, bValue, pmiValue
   return(aDictValue)


# Gain Ratio value  return a Dict (word:chi) from a corpusDFDict, subcorpusDFDict, 
def computeGR(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue  = float(corpusSize1) + corpusSize2
   acValue = float(corpusSize1)
   bdValue = float(corpusSize2)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0
      abValue = aValue + bValue
      cValue = acValue - aValue
      dValue = bdValue - bValue
      part1 = ((aValue) * nValue) / (abValue * acValue)
      part2 = (cValue * nValue) / (acValue * (cValue + dValue))
      if (aValue > 0):
         grValue =  (aValue / nValue) * math.log(part1,2)
      else:
         grValue = 0.0
      if (cValue > 0):
         grValue = grValue +  (cValue/nValue) * math.log(part2,2)
      aDictValue [aKey] = grValue
#      print aKey, aValue, bValue, cValue, dValue, grValue
   return(aDictValue)


# Information gain value  return a Dict (word:chi) from a corpusDFDict, subcorpusDFDict, 
def computeIG(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue  = float(corpusSize1) + corpusSize2
   acValue = float(corpusSize1)
   bdValue = float(corpusSize2)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0
      abValue = aValue + bValue
      cValue  = acValue - aValue
      dValue  = bdValue - bValue

      if (aValue > 0):
         part1 = aValue / nValue
         part1 = part1 * math.log ((aValue * nValue) / (abValue * acValue))
      else:
         part1 = 0
      if (bValue > 0):        
         part2 = bValue / nValue
         part2 = part2 * math.log ((bValue * nValue) / (abValue * bdValue))
      else:
         part2 = 0
      if (cValue > 0):
         part3 = cValue / nValue
         part3 = part3 * math.log ((cValue * nValue) / (acValue * (cValue+dValue)))
      else:
         part3 = 0
      if (dValue > 0):
         part4 = dValue / nValue
         part4 = part4 * math.log ((dValue * nValue) / ((bdValue) * (cValue+dValue)))
      else:
         part4 = 0
      igValue = part1 + part2 + part3 + part4         
      aDictValue [aKey] = igValue
#      print aKey, aValue, bValue, cValue, dValue, igValue
   return(aDictValue)


def computeCC(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue = float(corpusSize1) + corpusSize2
   acValue = float(corpusSize1)
   bdValue = float(corpusSize2)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = float(aDict2DF[aKey])
      else:
         bValue = 0.0
      if (aDict1DF.has_key(aKey)):
         aValue = float(aDict1DF[aKey])
      else:
         aValue = 0.0
      abValue = aValue + bValue
      cValue = acValue - aValue
      dValue = bdValue - bValue
      num = math.sqrt(nValue) * ((aValue * dValue)  - (cValue * bValue))
      denum = math.sqrt(abValue * acValue * bdValue * (cValue + dValue))
      if ((num <= 0) | (denum <= 0)):
         ccValue = 0.0
      else:
         ccValue =   num / denum
      aDictValue[aKey] = ccValue
#      print aKey, aValue, bValue, cValue, dValue, ccValue
   return(aDictValue)

# Relevance score  return a Dict 
def computeRS(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue = float(corpusSize1 + corpusSize2)
   acValue = float(corpusSize1)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = aDict2DF[aKey]
      else:
         bValue = 0
      if (aDict1DF.has_key(aKey)):
         aValue = aDict1DF[aKey]
      else:
         aValue = 0
      cValue = acValue - aValue
      dValue = nValue - acValue - bValue
      num = (float(aValue) / acValue) + 1.0
      denum = (float(bValue) / (bValue + dValue)) + 1.0
      rsValue = math.log ( (num / denum), 2) 
      aDictValue [aKey] = rsValue
#      print aKey, aValue, bValue, cValue, dValue, rsValue
   return(aDictValue)


# WLLR Eighted Log LikelihoodRatio Measure  return a Dict (word:wllr) from a corpusDFDict, subcorpusDFDict, 
def computeWLLR(aDict1DF, aDict2DF, corpusSize1, corpusSize2, aSet=[]):
   aDictValue = {}
   nValue = float(corpusSize1 + corpusSize2)
   acValue = float(corpusSize1)
   if (len(aSet) == 0):
      aKeySet = union(aDict1DF.keys(), aDict2DF.keys())
   else:
      aKeySet = aSet
   for aKey in aKeySet:
      if (aDict2DF.has_key(aKey)):
         bValue = aDict2DF[aKey]
      else:
         bValue = 0.0001
      if (aDict1DF.has_key(aKey)):
         aValue = aDict1DF[aKey]
      else:
         aValue = 0
      cValue = acValue - aValue
      dValue = nValue - acValue - bValue
      num = float(aValue) / acValue
      denum = float(bValue) / (bValue + dValue)
      if (num <= 0) :
         wllrValue = -1.0
      elif (denum <= 0):
         wllrValue = 1.0 + (float(aValue) / nValue) 
      else:
         wllrValue =  num *  math.log ( (num / denum), 2) 
      aDictValue[aKey] = wllrValue
#      print aKey, aValue, bValue, cValue, dValue, wllrValue
   return(aDictValue)



#
# Read files 
#

def createDictionary(aFileName):
   myFile = codecs.open(aFileName, encoding='utf-8', mode="r")
   aDict = {}
   aNumber = 0
   aLine = myFile.readline()
   myLine = string.split(aLine)
   myLine = myLine[3:]
   while aLine:
      aNumber = aNumber + 1
      for aWord in myLine:
         if (aDict.has_key(aWord)):
            aValue = aDict[aWord]
         else:
            aValue = 0
         aDict[aWord] = aValue + 1
      aLine = myFile.readline()
      myLine = string.split(aLine)
      myLine = myLine[3:]
   myFile.close()
   return(aDict, aNumber)

def createDFDictionary(aFileName, startPos=3):
   myFile = codecs.open(aFileName, encoding='utf-8', mode="r")
   aDict = {}
   aNumber = 0
   aLine = myFile.readline()
   myLine = string.split(aLine)
   myLine = myLine[startPos:]
   while aLine:
      aDict1 = {}
      aNumber = aNumber + 1
      for aWord in myLine:
         aDict1[aWord] = 1
      for aWord in aDict1.keys():
         if (aDict.has_key(aWord)):
            aValue = aDict[aWord]
         else:
            aValue = 0
         aDict[aWord] = aValue + 1
      aLine = myFile.readline()
      myLine = string.split(aLine)
      myLine = myLine[startPos:]
   myFile.close()
   return(aDict, aNumber)


def readDictionary(aFileName):
   myFile = codecs.open(aFileName, encoding='utf-8', mode="r")
   aDict = {}
   aLine = myFile.readline()
   aLine = aLine.replace(u"\n", "")
   aPos = aLine.find(u':')
   aKey = aLine[:aPos]
   aValue = aLine[aPos+1:]
   while aLine:
      aDict[aValue] = float(aKey)
      aLine = myFile.readline()
      aLine = aLine.replace(u"\n", "")
      aPos = aLine.find(u':')
      aKey = aLine[:aPos]
      aValue = aLine[aPos+1:]
   myFile.close()
   return(aDict)




#
# Union and intersection
#

def intersection (*args):
   res = []
   for x in args[0]:
      for other in args[1:]:
         if x not in other: break
      else:
        res.append(x)
   return res

def union (*args):
   res = []
   for seq in args:
      for x in seq:
         if x not in res:
            res.append(x)
   return res

#
# Manipulating a dictionary
#

# remove entries for which the tf value <= minTF
def filterDictTF(aDict, minTF):
   for aKey in aDict.keys():
      if (aDict[aKey] <= minTF):
         del aDict[aKey]

# generate a new dict by normalizing the values (tf -> ntf)
# use to create a dict with relative frequencies
def reverseDict(aDict):
   aNewDict = {}
   for aKey in aDict.keys():
      aValue = aDict[aKey]
      if (aNewDict.has_key(aValue)):
         aList = aNewDict[aValue]
      else:
         aList = []
      aList.append(aKey)
      aNewDict[aValue] = aList
   return(aNewDict)

def dumpDictionary(aDictA, aFileName=""):
   aDict = reverseDict(aDictA)
   if (len(aFileName) > 0):
      myFile = codecs.open(aFileName, encoding='utf-8', mode="a")
   for aKey in sorted(aDict.keys(),None, None, True):
      aListKeyA = aDict[aKey]
      for aKeyA in aListKeyA:
         aValueA = aDictA[aKeyA]
         if (len(aFileName) > 0):
            myFile.write(str(aValueA)+":"+aKeyA+"\n")
         else:
            print str(aValueA)+"  :  "+aKeyA
   if (len(aFileName) > 0):
      myFile.close()
   return(aDict)

def extractTopSortedDict(aDict, aLimit, positiveOnly=False):
   if not(positiveOnly):
      return(extractTopSortedDictPrivate(aDict, aLimit))
   aNewDict = reverseDict(aDict)
   printed = 0
   aTopDict = {}
   for aKey in sorted(aNewDict.keys(),None,None,True):
      if (printed < aLimit):
         for aWord in aNewDict[aKey]:
            if (aKey < 0):
               break
            aTopDict[aWord] = aKey
            printed = printed + 1
            if (printed > aLimit):
               break
      else:
         break
   return(aTopDict)

def extractTopSortedDictPrivate(aDict, aLimit):
   aNewDict = reverseDict(aDict)
   printed = 0
   aTopDict = {}
   for aKey in sorted(aNewDict.keys(),None,None,True):
      if (printed < aLimit):
         for aWord in aNewDict[aKey]:
            aTopDict[aWord] = aKey
            printed = printed + 1
            if (printed > aLimit):
               break
      else:
         break
   return(aTopDict)


def extractBottomSortedDict(aDict, aLimit):
   aNewDict = reverseDict(aDict)
   printed = 0
   aBottomDict = {}
   for aKey in sorted(aNewDict.keys()):
      if (printed < aLimit):
         for aWord in aNewDict[aKey]:
            aBottomDict[aWord] = aKey
            printed = printed + 1
            if (printed > aLimit):
               break
      else:
         break
   return(aBottomDict)

# print only when not imported
if __name__ == '__main__':
   print "The module TermSelection.py was uploaded\n"
