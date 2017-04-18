from SimPEG import Mesh, Regularization, Maps, Utils, EM
from SimPEG.EM.Static import DC
import numpy as np
import matplotlib.pyplot as plt 
#%matplotlib inline
import copy
#import pandas as pd
#from scipy.sparse import csr_matrix, spdiags, dia_matrix,diags
#from scipy.sparse.linalg import spsolve
from  scipy.stats import norm,multivariate_normal
import sys
path ="../pymatsolver/" 
path = "../../../Documents/pymatsolver/"
sys.path.append(path)
from pymatsolver import PardisoSolver
#from scipy.interpolate import LinearNDInterpolator, interp1d
#from sklearn.mixture import GaussianMixture
from SimPEG import DataMisfit, Regularization, Optimization, InvProblem, Directives, Inversion
import SimPEG
import scipy.sparse as sp

import os
import glob

#Remove older results
files = glob.glob('./*.npz')
for f in files:
    os.remove(f)

#2D model
csx, csy, csz = 0.25,0.25,0.25
# Number of core cells in each directiPon s
ncx, ncz = 123,41
# Number of padding cells to add in each direction
npad = 12
# Vectors of cell lengthts in each direction
hx = [(csx,npad, -1.5),(csx,ncx),(csx,npad, 1.5)]
hz= [(csz,npad,-1.5),(csz,ncz)]
# Create mesh
mesh = Mesh.TensorMesh([hx, hz],x0="CN")
# Map mesh coordinates from local to UTM coordiantes
#mesh.x0[2] = mesh.x0[2]-mesh.vectorCCz[-npad-1]
mesh.x0[1] = mesh.x0[1]+csz/2.
#mesh.x0[0] = mesh.x0[0]+csx/2.
#mesh.plotImage(np.ones(mesh.nC)*np.nan, grid=True)
#mesh.plotImage(np.ones(mesh.nC)*np.nan, grid=True)
#plt.gca().set_xlim([-20,20])
#plt.gca().set_ylim([-15,0])
#mesh.plotGrid()
#plt.gca().set_aspect('equal')
#plt.show()

print "Mesh Size: ", mesh.nC


#Model Creation
lnsig_air = 1e-8;

x0,z0, r0 = -6., -4., 3.
x1,z1, r1 =  6., -4., 3.

ln_sigback = -5.
ln_sigc = -3.
ln_sigr = -7.

noisemean = 0.
noisevar = 0.0
overburden_extent = 0.
ln_over = -4.
#m = (lnsig_background)*np.ones(mesh.nC);
#mu =np.ones(mesh.nC);
mtrue = ln_sigback*np.ones(mesh.nC) + norm(noisemean,noisevar).rvs(mesh.nC)

overb = (mesh.gridCC[:,1] >-overburden_extent) & (mesh.gridCC[:,1]<=0)
mtrue[overb] = ln_over*np.ones_like(mtrue[overb])+ norm(noisemean,noisevar).rvs(np.prod((mtrue[overb]).shape))


csph = (np.sqrt((mesh.gridCC[:,1]-z0)**2.+(mesh.gridCC[:,0]-x0)**2.))< r0
mtrue[csph] = ln_sigc*np.ones_like(mtrue[csph]) + norm(noisemean,noisevar).rvs(np.prod((mtrue[csph]).shape))

#Define the sphere limit
rsph = (np.sqrt((mesh.gridCC[:,1]-z1)**2.+(mesh.gridCC[:,0]-x1)**2.))< r1
mtrue[rsph] = ln_sigr*np.ones_like(mtrue[rsph]) + norm(noisemean,noisevar).rvs(np.prod((mtrue[rsph]).shape))

mtrue = Utils.mkvc(mtrue);

mesh.plotGrid()
plt.gca().set_xlim([-10,10])
plt.gca().set_ylim([-10,0])
xyzlim = np.r_[[[-10.,10.],[-10.,1.]]]
actind, meshCore = Utils.meshutils.ExtractCoreMesh(xyzlim,mesh)

plt.hist(mtrue[actind],bins =50,normed=True);

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
mm = meshCore.plotImage(mtrue[actind],ax = ax0)
plt.colorbar(mm[0])
ax0.set_aspect("equal")
#plt.show()

#Gradient array 1 2D
srclist = []
nSrc = 23
lines = 1
ylines = np.r_[0.]
xlines = np.r_[0.]
z = 0.

#xline
for k in range(lines):
    for i in range(nSrc):
        if i<=11:
            locA = np.r_[-14.+1., z]
            locB = np.r_[-8.+2.*i-1., z]
            #M = np.c_[np.arange(-12.,-12+2*(i+1),2),np.ones(i+1)*z]
            #N = np.c_[np.arange(-10.,-10+2*(i+1),2),np.ones(i+1)*z]
            M = np.c_[np.arange(-12.,10+1,2),np.ones(12)*z]
            N = np.c_[np.arange(-10.,12+1,2),np.ones(12)*z]
            rx = DC.Rx.Dipole(M,N)
            src= DC.Src.Dipole([rx],locA,locB)
            srclist.append(src)
            #print locA,locB,"\n",[M,N],"\n"
            
            #rx = DC.Rx.Dipole(-M,-N)
            #src= DC.Src.Dipole([rx],-locA,-locB)
            #srclist.append(src)
            #print -locA,-locB,"\n",[-M,-N],"\n"


            
        else:
            locA = np.r_[-14.+2*(i-11)+1., z]
            locB = np.r_[14.-1.,z]
            #M = np.c_[np.arange(locA[0]+1.,12.,2),np.ones(nSrc-i)*z]
            #N = np.c_[np.arange(locA[0]+3.,14.,2),np.ones(nSrc-i)*z]
            M = np.c_[np.arange(-12.,10+1,2),np.ones(12)*z]
            N = np.c_[np.arange(-10.,12+1,2),np.ones(12)*z]
            rx = DC.Rx.Dipole(M,N)
            src= DC.Src.Dipole([rx],locA,locB)
            srclist.append(src)
            #print "line2",locA,locB,"\n",[M,N],"\n"
            
            #rx = DC.Rx.Dipole(-M,-N)
            #src= DC.Src.Dipole([rx],-locA,-locB)
            #srclist.append(src)

mapping = Maps.ExpMap(mesh)
survey = DC.Survey(srclist)
problem = DC.Problem3D_CC(mesh, sigmaMap=mapping)
problem.pair(survey)
problem.Solver = PardisoSolver
survey.dobs = survey.dpred(mtrue)
survey.std = 0.05*np.ones_like(survey.dobs)
survey.eps = 1e-5*np.linalg.norm(survey.dobs)
print '# of data: ', survey.dobs.shape

class SimultaneousSrc(DC.Src.BaseSrc):
    """
    Dipole source
    """
    QW = None
    Q = None
    W = None
    def __init__(self, rxList,Q,W, **kwargs):
        
        SimPEG.Survey.BaseSrc.__init__(self, rxList, **kwargs)

    def eval(self, prob):
        return self.QW

class SimultaneousRx(DC.Rx.BaseRx):
    """
    SimultaneousRx receiver
    """

    def __init__(self, locs, rxType='phi', **kwargs):
        # We may not need this ...
        SimPEG.Survey.BaseRx.__init__(self, locs, rxType)
                    
    @property
    def nD(self):
        """Number of data in the receiver."""
        return self.locs.shape[0]

        # Not sure why ...
        # return int(self.locs[0].size / 2)

    def getP(self, mesh, Gloc):
        return self.locs

P = []
M = np.c_[np.arange(-12.,10+1,2),np.ones(12)*z]
N = np.c_[np.arange(-10.,12+1,2),np.ones(12)*z]
rx = DC.Rx.Dipole(M,N)
P = rx.getP(mesh,'CC')


#Update W Inversion
nsubSrc = 5
m0 = (-5.)*np.ones(mapping.nP);
miter = m0
n_its = 50
InnerIt = 3
dmisfitsub = []
dmisfitall = []
#beta schedule
beta = 1.
betalist = [beta]
coolingFactor = 2.
coolingRate = 3
W = np.random.randn(survey.nSrc,nsubSrc)

dmisAll = DataMisfit.l2_DataMisfit(survey)
dmisfitall.append(dmisAll.eval(m0)/survey.nD)
print "Starting Model Dmisfit compared to full dataset: ",dmisAll.eval(m0)/survey.nD
print "Check misfit with true model: ",dmisAll.eval(mtrue)/survey.nD

for it in range(n_its): 
    problem.unpair()
    problem.pair(survey)
    Q = problem.getRHS()
    sub = problem.getRHS().dot(W)

    rx_r = SimultaneousRx(locs=P)
    srcList_r = []
    for isrc in range(sub.shape[1]):
        src_r = SimultaneousSrc([rx_r], Q=Q[:,isrc],W=W[:,isrc],QW =Q.dot(W)[:,isrc])
        srcList_r.append(src_r)
    survey_r = DC.Survey(srcList_r)

    problem.unpair()
    problem.pair(survey_r)

    d = survey_r.dpred(mtrue)
    survey_r.dobs = d
    survey_r.std = np.ones_like(d)*0.05
    survey_r.eps = 1e-5*np.linalg.norm(survey_r.dobs)

    print '# of data: ', survey_r.dobs.shape

    regmesh = mesh;
    dmis = DataMisfit.l2_DataMisfit(survey_r)
    reg = Regularization.Tikhonov(regmesh)#,mapping = mapping)#,indActive=actind)
    reg.mref = m0
    opt = Optimization.InexactGaussNewton(maxIter=1,tolX=1e-6)
    opt.remember('xc')
    invProb = InvProblem.BaseInvProblem(dmis, reg, opt)
    #beta = Directives.BetaEstimate_ByEig(beta0= 10.,beta0_ratio=1e0)
    reg.alpha_s = 1e-6;
    invProb.beta = beta
    #betaSched = Directives.BetaSchedule(coolingFactor=5, coolingRate=2)
    #sav0 = Directives.SaveEveryIteration()
    #sav1 = Directives.SaveModelEveryIteration()
    #sav2 = Directives.SaveOutputDictEveryIteration()
    inv = Inversion.BaseInversion(invProb)#, directiveList=[sav2])#[beta,betaSched])#sav0,sav1,

    msimple = inv.run(miter);
    beta = invProb.beta
    if np.mod(it+1,coolingRate) ==0:
        beta = beta/coolingFactor
    betalist.append(beta)
    miter = copy.deepcopy(msimple)
    dmisfitsub.append(dmis.eval(msimple)/survey_r.nD)
    print "Dmisfit compared to sub dataset: ",dmis.eval(msimple)/survey_r.nD
    print "Check misfit with true model: ",dmis.eval(mtrue)/survey_r.nD

    problem.unpair()
    problem.pair(survey)
    dmisAll = DataMisfit.l2_DataMisfit(survey)
    dmisfitall.append(dmisAll.eval(msimple)/survey.nD)
    print "Dmisfit compared to full dataset: ",dmisAll.eval(msimple)/survey.nD
    print "Check misfit with true model: ",dmisAll.eval(mtrue)/survey.nD


    if np.mod(it+1,InnerIt) ==0:
        W = np.random.randn(survey.nSrc,nsubSrc)
        print 'update W'
    

    #mm = mesh.plotImage(miter)
    #plt.colorbar(mm[0])
    #plt.gca().set_xlim([-10.,10.])
    #plt.gca().set_ylim([-10.,0.])

np.save('./dmisfitsub.npy',dmisfitsub)
np.save('./dmisfitall.npy',dmisfitall)
np.save('./beta.npy',betalist)
np.save('./finalresult',msimple)

#plt.show()'

