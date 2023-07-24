import ROOT
import argparse
import os
import math
from array import array

def correctedMET (originalMet, originalMet_phi, npv, runnb, isMC, yeara, isUL =True, ispuppi=False) :

  year=str(yeara)
  #print 'some info=================================', isMC, year, isUL, ispuppi
  if(npv>100) :npv=100
  runera =-1
  usemetv2 =False
  year = year.replace("postVFP", "nonAPV")
  year = year.replace("preVFP", "APV")

  if (isMC and year == "2016" and isUL) :  runera ="yUL2016MCnonAPV"
  elif(isMC and year == "2016APV" and isUL): runera = "yUL2016MCAPV"
  elif(isMC and year == "2016nonAPV" and isUL): runera = "yUL2016MCnonAPV"
  elif(isMC and year == "2017" and isUL) :runera = "yUL2017MC"
  elif(isMC and year == "2018" and isUL) :runera = "yUL2018MC"
  elif(not isMC and runnb >=315252 and runnb <=316995 and isUL): runera = "yUL2018A"
  elif(not isMC and runnb >=316998 and runnb <=319312 and isUL): runera = "yUL2018B"
  elif(not isMC and runnb >=319313 and runnb <=320393 and isUL): runera = "yUL2018C"
  elif(not isMC and runnb >=320394 and runnb <=325273 and isUL): runera = "yUL2018D"

  elif(not isMC and runnb >=297020 and runnb <=299329 and isUL):
      runera = "yUL2017B"
      usemetv2 =False
  elif(not isMC and runnb >=299337 and runnb <=302029 and isUL):
      runera = "yUL2017C"
      usemetv2 =False
  elif(not isMC and runnb >=302030 and runnb <=303434 and isUL):
      runera = "yUL2017D"
      usemetv2 =False
  elif(not isMC and runnb >=303435 and runnb <=304826 and isUL):
      runera = "yUL2017E"
      usemetv2 =False
  elif(not isMC and runnb >=304911 and runnb <=306462 and isUL):
      runera = "yUL2017F"
      usemetv2 =False

  elif(not isMC and runnb >=272007 and runnb <=275376 and isUL): runera = "yUL2016B"
  elif(not isMC and runnb >=275657 and runnb <=276283 and isUL): runera = "yUL2016C"
  elif(not isMC and runnb >=276315 and runnb <=276811 and isUL): runera = "yUL2016D"
  elif(not isMC and runnb >=276831 and runnb <=277420 and isUL): runera = "yUL2016E"
  elif(not isMC and ((runnb >=277772 and runnb <=278768) or runnb==278770) and isUL): runera = "yUL2016F"
  elif(not isMC and ((runnb >=278801 and runnb <=278808) or runnb==278769) and isUL): runera = "yUL2016Flate"
  elif(not isMC and runnb >=278820 and runnb <=280385 and isUL): runera = "yUL2016G"
  elif(not isMC and runnb >=280919 and runnb <=284044 and isUL): runera = "yUL2016H"

  else :   
      print 'failed, ============================================>', year, runera, isMC
      return [1, 1, originalMet, originalMet_phi]
  #print 'info', runera, isMC
 
  #print '============================================>', year, runera, isMC
  
  ''' 
  if(isMC and year == "2016" and !isUL) :runera = "y2016MC"
  elif(isMC and year == "2017" and !isUL) : 
      runera = "y2017MC"
      usemetv2 =True
  elif(isMC and year == "2018" and !isUL): runera = "y2018MC"
  elif(not isMC and runnb >=272007 and runnb <=275376 and !isUL): runera = "y2016B"
  elif(not isMC and runnb >=275657 and runnb <=276283 and !isUL): runera = "y2016C"
  elif(not isMC and runnb >=276315 and runnb <=276811 and !isUL): runera = "y2016D"
  elif(not isMC and runnb >=276831 and runnb <=277420 and !isUL): runera = "y2016E"
  elif(not isMC and runnb >=277772 and runnb <=278808 and !isUL): runera = "y2016F"
  elif(not isMC and runnb >=278820 and runnb <=280385 and !isUL): runera = "y2016G"
  elif(not isMC and runnb >=280919 and runnb <=284044 and !isUL): runera = "y2016H"
  
  elif(not isMC and runnb >=297020 and runnb <=299329 and !isUL): 
      runera = "y2017B"
      usemetv2 =True
  elif(not isMC and runnb >=299337 and runnb <=302029 and !isUL):
      runera = "y2017C"
      usemetv2 =True
  elif(not isMC and runnb >=302030 and runnb <=303434 and !isUL):
      runera = "y2017D"
      usemetv2 =True
  elif(not isMC and runnb >=303435 and runnb <=304826 and !isUL):
      runera = "y2017E"
      usemetv2 =True
  elif(not isMC and runnb >=304911 and runnb <=306462 and !isUL):
      runera = "y2017F"
      usemetv2 =True
  elif(not isMC and runnb >=315252 and runnb <=316995 and !isUL): runera = "y2018A"
  elif(not isMC and runnb >=316998 and runnb <=319312 and !isUL): runera = "y2018B"
  elif(not isMC and runnb >=319313 and runnb <=320393 and !isUL): runera = "y2018C"
  elif(not isMC and runnb >=320394 and runnb <=325273 and !isUL): runera = "y2018D"
  '''

  
  
  METxcorr,METycorr=0.,0.

  if not usemetv2:
    if not ispuppi:
	if(runera=="y2016B"): METxcorr = -(-0.0478335*npv -0.108032)
	if(runera=="y2016B"): METycorr = -(0.125148*npv +0.355672)
	if(runera=="y2016C"): METxcorr = -(-0.0916985*npv +0.393247)
	if(runera=="y2016C"): METycorr = -(0.151445*npv +0.114491)
	if(runera=="y2016D"): METxcorr = -(-0.0581169*npv +0.567316)
	if(runera=="y2016D"): METycorr = -(0.147549*npv +0.403088)
	if(runera=="y2016E"): METxcorr = -(-0.065622*npv +0.536856)
	if(runera=="y2016E"): METycorr = -(0.188532*npv +0.495346)
	if(runera=="y2016F"): METxcorr = -(-0.0313322*npv +0.39866)
	if(runera=="y2016F"): METycorr = -(0.16081*npv +0.960177)
	if(runera=="y2016G"): METxcorr = -(0.040803*npv -0.290384)
	if(runera=="y2016G"): METycorr = -(0.0961935*npv +0.666096)
	if(runera=="y2016H"): METxcorr = -(0.0330868*npv -0.209534)
	if(runera=="y2016H"): METycorr = -(0.141513*npv +0.816732)
	if(runera=="y2017B"): METxcorr = -(-0.259456*npv +1.95372)
	if(runera=="y2017B"): METycorr = -(0.353928*npv -2.46685)
	if(runera=="y2017C"): METxcorr = -(-0.232763*npv +1.08318)
	if(runera=="y2017C"): METycorr = -(0.257719*npv -1.1745)
	if(runera=="y2017D"): METxcorr = -(-0.238067*npv +1.80541)
	if(runera=="y2017D"): METycorr = -(0.235989*npv -1.44354)
	if(runera=="y2017E"): METxcorr = -(-0.212352*npv +1.851)
	if(runera=="y2017E"): METycorr = -(0.157759*npv -0.478139)
	if(runera=="y2017F"): METxcorr = -(-0.232733*npv +2.24134)
	if(runera=="y2017F"): METycorr = -(0.213341*npv +0.684588)
	if(runera=="y2018A"): METxcorr = -(0.362865*npv -1.94505)
	if(runera=="y2018A"): METycorr = -(0.0709085*npv -0.307365)
	if(runera=="y2018B"): METxcorr = -(0.492083*npv -2.93552)
	if(runera=="y2018B"): METycorr = -(0.17874*npv -0.786844)
	if(runera=="y2018C"): METxcorr = -(0.521349*npv -1.44544)
	if(runera=="y2018C"): METycorr = -(0.118956*npv -1.96434)
	if(runera=="y2018D"): METxcorr = -(0.531151*npv -1.37568)
	if(runera=="y2018D"): METycorr = -(0.0884639*npv -1.57089)
	if(runera=="y2016MC"): METxcorr = -(-0.195191*npv -0.170948)
	if(runera=="y2016MC"): METycorr = -(-0.0311891*npv +0.787627)
	if(runera=="y2017MC"): METxcorr = -(-0.217714*npv +0.493361)
	if(runera=="y2017MC"): METycorr = -(0.177058*npv -0.336648)
	if(runera=="y2018MC"): METxcorr = -(0.296713*npv -0.141506)
	if(runera=="y2018MC"): METycorr = -(0.115685*npv +0.0128193)

	#//UL2017
	if(runera=="yUL2017B"): METxcorr = -(-0.211161*npv +0.419333)
	if(runera=="yUL2017B"): METycorr = -(0.251789*npv +-1.28089)
	if(runera=="yUL2017C"): METxcorr = -(-0.185184*npv +-0.164009)
	if(runera=="yUL2017C"): METycorr = -(0.200941*npv +-0.56853)
	if(runera=="yUL2017D"): METxcorr = -(-0.201606*npv +0.426502)
	if(runera=="yUL2017D"): METycorr = -(0.188208*npv +-0.58313)
	if(runera=="yUL2017E"): METxcorr = -(-0.162472*npv +0.176329)
	if(runera=="yUL2017E"): METycorr = -(0.138076*npv +-0.250239)
	if(runera=="yUL2017F"): METxcorr = -(-0.210639*npv +0.72934)
	if(runera=="yUL2017F"): METycorr = -(0.198626*npv +1.028)
	if(runera=="yUL2017MC"): METxcorr = -(-0.300155*npv +1.90608)
	if(runera=="yUL2017MC"): METycorr = -(0.300213*npv +-2.02232)

	#//UL2018
	if(runera=="yUL2018A"): METxcorr = -(0.263733*npv +-1.91115)
	if(runera=="yUL2018A"): METycorr = -(0.0431304*npv +-0.112043)
	if(runera=="yUL2018B"): METxcorr = -(0.400466*npv +-3.05914)
	if(runera=="yUL2018B"): METycorr = -(0.146125*npv +-0.533233)
	if(runera=="yUL2018C"): METxcorr = -(0.430911*npv +-1.42865)
	if(runera=="yUL2018C"): METycorr = -(0.0620083*npv +-1.46021)
	if(runera=="yUL2018D"): METxcorr = -(0.457327*npv +-1.56856)
	if(runera=="yUL2018D"): METycorr = -(0.0684071*npv +-0.928372)
	if(runera=="yUL2018MC"): METxcorr = -(0.183518*npv +0.546754)
	if(runera=="yUL2018MC"): METycorr = -(0.192263*npv +-0.42121)

	#//UL2016
	if(runera=="yUL2016B"): METxcorr = -(-0.0214894*npv +-0.188255)
	if(runera=="yUL2016B"): METycorr = -(0.0876624*npv +0.812885)
	if(runera=="yUL2016C"): METxcorr = -(-0.032209*npv +0.067288)
	if(runera=="yUL2016C"): METycorr = -(0.113917*npv +0.743906)
	if(runera=="yUL2016D"): METxcorr = -(-0.0293663*npv +0.21106)
	if(runera=="yUL2016D"): METycorr = -(0.11331*npv +0.815787)
	if(runera=="yUL2016E"): METxcorr = -(-0.0132046*npv +0.20073)
	if(runera=="yUL2016E"): METycorr = -(0.134809*npv +0.679068)
	if(runera=="yUL2016F"): METxcorr = -(-0.0543566*npv +0.816597)
	if(runera=="yUL2016F"): METycorr = -(0.114225*npv +1.17266)
	if(runera=="yUL2016Flate"): METxcorr = -(0.134616*npv +-0.89965)
	if(runera=="yUL2016Flate"): METycorr = -(0.0397736*npv +1.0385)
	if(runera=="yUL2016G"): METxcorr = -(0.121809*npv +-0.584893)
	if(runera=="yUL2016G"): METycorr = -(0.0558974*npv +0.891234)
	if(runera=="yUL2016H"): METxcorr = -(0.0868828*npv +-0.703489)
	if(runera=="yUL2016H"): METycorr = -(0.0888774*npv +0.902632)
	if(runera=="yUL2016MCnonAPV"): METxcorr = -(-0.153497*npv +-0.231751)
	if(runera=="yUL2016MCnonAPV"): METycorr = -(0.00731978*npv +0.243323)
	if(runera=="yUL2016MCAPV"): METxcorr = -(-0.188743*npv +0.136539)
	if(runera=="yUL2016MCAPV"): METycorr = -(0.0127927*npv +0.117747)


    if ispuppi :

	if(runera=="yUL2017B"): METxcorr = -(-0.00382117*npv +-0.666228)
	if(runera=="yUL2017B"): METycorr = -(0.0109034*npv +0.172188)
	if(runera=="yUL2017C"): METxcorr = -(-0.00110699*npv +-0.747643)
	if(runera=="yUL2017C"): METycorr = -(-0.0012184*npv +0.303817)
	if(runera=="yUL2017D"): METxcorr = -(-0.00141442*npv +-0.721382)
	if(runera=="yUL2017D"): METycorr = -(-0.0011873*npv +0.21646)
	if(runera=="yUL2017E"): METxcorr = -(0.00593859*npv +-0.851999)
	if(runera=="yUL2017E"): METycorr = -(-0.00754254*npv +0.245956)
	if(runera=="yUL2017F"): METxcorr = -(0.00765682*npv +-0.945001)
	if(runera=="yUL2017F"): METycorr = -(-0.0154974*npv +0.804176)
	if(runera=="yUL2017MC"): METxcorr = -(-0.0102265*npv +-0.446416)
	if(runera=="yUL2017MC"): METycorr = -(0.0198663*npv +0.243182)

	#//UL2018Puppi
	if(runera=="yUL2018A"): METxcorr = -(-0.0073377*npv +0.0250294)
	if(runera=="yUL2018A"): METycorr = -(-0.000406059*npv +0.0417346)
	if(runera=="yUL2018B"): METxcorr = -(0.00434261*npv +0.00892927)
	if(runera=="yUL2018B"): METycorr = -(0.00234695*npv +0.20381)
	if(runera=="yUL2018C"): METxcorr = -(0.00198311*npv +0.37026)
	if(runera=="yUL2018C"): METycorr = -(-0.016127*npv +0.402029)
	if(runera=="yUL2018D"): METxcorr = -(0.00220647*npv +0.378141)
	if(runera=="yUL2018D"): METycorr = -(-0.0160244*npv +0.471053)
	if(runera=="yUL2018MC"): METxcorr = -(-0.0214557*npv +0.969428)
	if(runera=="yUL2018MC"): METycorr = -(0.0167134*npv +0.199296)

	#//UL2016Puppi
	if(runera=="yUL2016B"): METxcorr = -(-0.00109025*npv +-0.338093)
	if(runera=="yUL2016B"): METycorr = -(-0.00356058*npv +0.128407)
	if(runera=="yUL2016C"): METxcorr = -(-0.00271913*npv +-0.342268)
	if(runera=="yUL2016C"): METycorr = -(0.00187386*npv +0.104)
	if(runera=="yUL2016D"): METxcorr = -(-0.00254194*npv +-0.305264)
	if(runera=="yUL2016D"): METycorr = -(-0.00177408*npv +0.164639)
	if(runera=="yUL2016E"): METxcorr = -(-0.00358835*npv +-0.225435)
	if(runera=="yUL2016E"): METycorr = -(-0.000444268*npv +0.180479)
	if(runera=="yUL2016F"): METxcorr = -(0.0056759*npv +-0.454101)
	if(runera=="yUL2016F"): METycorr = -(-0.00962707*npv +0.35731)
	if(runera=="yUL2016Flate"): METxcorr = -(0.0234421*npv +-0.371298)
	if(runera=="yUL2016Flate"): METycorr = -(-0.00997438*npv +0.0809178)
	if(runera=="yUL2016G"): METxcorr = -(0.0182134*npv +-0.335786)
	if(runera=="yUL2016G"): METycorr = -(-0.0063338*npv +0.093349)
	if(runera=="yUL2016H"): METxcorr = -(0.015702*npv +-0.340832)
	if(runera=="yUL2016H"): METycorr = -(-0.00544957*npv +0.199093)
	if(runera=="yUL2016MCnonAPV"): METxcorr = -(-0.0058341*npv +-0.395049)
	if(runera=="yUL2016MCnonAPV"): METycorr = -(0.00971595*npv +-0.101288)
	if(runera=="yUL2016MCAPV"): METxcorr = -(-0.0060447*npv +-0.4183)
	if(runera=="yUL2016MCAPV"): METycorr = -(0.008331*npv +-0.0990046)



	

  #print "METcorr: ",runnb,runera,npv,METxcorr,METycorr

  CorrectedMET_x = originalMet * math.cos(originalMet_phi) + METxcorr
  CorrectedMET_y = originalMet * math.sin(originalMet_phi) + METycorr

  CorrectedMET = math.sqrt(CorrectedMET_x*CorrectedMET_x + CorrectedMET_y*CorrectedMET_y)
  CorrectedMETPhi = 0.

  if CorrectedMET_x==0 and CorrectedMET_y>0 : CorrectedMETPhi = math.pi
  elif CorrectedMET_x==0 and CorrectedMET_y<0  :CorrectedMETPhi = -math.pi
  elif CorrectedMET_x >0 : CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x)
  elif CorrectedMET_x <0 and CorrectedMET_y>0 : CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x) + math.pi
  elif CorrectedMET_x <0 and CorrectedMET_y<0 : CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x) - math.pi  
  else: CorrectedMETPhi =0.

  return [CorrectedMET_x, CorrectedMET_y, CorrectedMET, CorrectedMETPhi]



