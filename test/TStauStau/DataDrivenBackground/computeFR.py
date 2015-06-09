
def extrapolateFR(rq1, rg1, f1, rq2, rg2, f2, rq, rg, verbose = 0):
  if(verbose > 1):
    print " In the first  region: rq=", rq1, ", rg=", rg1, ", f=", f1 
    print " In the second region: rq=", rq2, ", rg=", rg2, ", f=", f2
    print " Extrapolating f in the region where: rq=", rq, ", rg=", rg 

  fg = (f2*rq1*(rq2+rg2) - f1*rq2*(rq1+rg1))/(rg2*rq1 - rg1*rq2)
  fq = (f1*rg2*(rq1+rg1) - f2*rg1*(rq2+rg2))/(rg2*rq1 - rg1*rq2)
  if(verbose > 0):
    print "  fq = ", fq
    print "  fg = ", fg

  fr = (fq*rq + fg*rg)/(rq+rg)
  if(verbose > 0):
    print "  f  = ", fr

  return fr

## Opposite Sign
rqOSWJets = 0.908356# +- 0.00442847
rgOSWjets = 0.0868501# +- 0.00102724
rqOStotal = 0.868817# +- 0.0279279
rgOStotal = 0.0974855# +- 0.0036059

fOSWJets = 0.507862# +- 0.00296636
fOSdata  = 0.584205# +- 0.00119172
fOSdatap = 0.51368# +- 0.00189922

print "____________________"
print "WJets OS test"
fq = 0.514
fg = 0.453
print "f(compu) = ", (fq*rqOSWJets + fg*rgOSWjets)/(rqOSWJets + rgOSWjets)
print "f(WJets) = ", fOSWJets

print "____________________"

## Same Sign
rqSSWJets = 0.802126# +- 0.0081367
rgSSWJets = 0.180924# +- 0.00315348
rqSStotal = 0.845346# +- 0.0749795
rgSStotal = 0.142974# +- 0.0142453

fSSWJets = 0.415558# +- 0.00539257
fSSdata  = 0.465788# +- 0.00169674
fSSdatap = 0.464945# +- 0.00172752

## Charge Symmetric
rqCSWJets = 0.876064# +- 0.0684528
rgCSWJets = 0.113009# +- 0.0167653

fCSWJets = 0.487688# +- 0.046073
fCSdata  = 0.312057# +- 0.0605102
fCSdatap = 0.46329# +- 0.00257595

## Inverted MET cut (MET < 30)
rqIMWJets = 0.910758# +- 0.00748617
rgIMWJets = 0.0850577# +- 0.00174246
rqIMtotal = 0.844293# +- 0.05715
rgIMtotal = 0.100065# +- 0.00675946

fIMWJets = 0.516165# +- 0.00523023
fIMdata  = 0.606565# +- 0.00166898
fIMdatap = 0.515584# +- 0.00312075

## MET Cut (MET > 30)  aka Analysis selection
rqASWJets = 0.907016# +- 0.00549148
rgASWJets = 0.0878495# +- 0.00127125
rqAStotal = 0.891147# +- 0.00875507
rgAStotal = 0.0951272# +- 0.00322537

fASWJets = 0.504593# +- 0.00366414
fASdata  = 0.558289# +- 0.00169775
fASdatap = 0.511737# +- 0.0022508



print "Using SS region and InvMet region to calculate FR, in WJets with WJets ratio"
extrapolateFR(rqSSWJets, rgSSWJets, fSSWJets, rqIMWJets, rgIMWJets, fIMWJets, rqASWJets, rgASWJets, verbose = 2)

print "Using SS region and InvMet region to calculate FR, in WJets with total ratio"
extrapolateFR(rqSStotal, rgSStotal, fSSWJets, rqIMtotal, rgIMtotal, fIMWJets, rqAStotal, rgAStotal, verbose = 2)

print "Using SS region and InvMet region to calculate FR, in data (subtracting prompt) with total ratio"
extrapolateFR(rqSStotal, rgSStotal, fSSdatap, rqIMtotal, rgIMtotal, fIMdatap, rqAStotal, rgAStotal, verbose = 2)


print "_____________________________________________________________________________________________________________"
print "Using SS region and CS region to calculate FR, in WJets with WJets ratio"
extrapolateFR(rqSSWJets, rgSSWJets, fSSWJets, rqCSWJets, rgCSWJets, fCSWJets, rqASWJets, rgASWJets, verbose = 2)

#print "Using SS region and CS region to calculate FR, in WJets with total ratio"
#extrapolateFR(rqSStotal, rgSStotal, fSSWJets, rqCStotal, rgCStotal, fCSWJets, rqAStotal, rgAStotal, verbose = 2)

#print "Using SS region and CS region to calculate FR, in data (subtracting prompt) with total ratio"
#extrapolateFR(rqSStotal, rgSStotal, fSSdatap, rqCStotal, rgCStotal, fCSdatap, rqAStotal, rgAStotal, verbose = 2)
