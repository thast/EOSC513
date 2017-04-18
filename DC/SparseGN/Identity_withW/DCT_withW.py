from SimPEG import Mesh, Regularization, Maps, Utils, EM
from SimPEG.EM.Static import DC
import numpy as np
import matplotlib.pyplot as plt 
#%matplotlib inline
import copy
import pandas as pd
from scipy.sparse import csr_matrix, spdiags, dia_matrix,diags
from scipy.sparse.linalg import spsolve
from  scipy.stats import norm,multivariate_normal
import sys
path ="../pymatsolver/" 
path = "../../../Documents/pymatsolver/"
sys.path.append(path)
from pymatsolver import PardisoSolver
from scipy.interpolate import LinearNDInterpolator, interp1d
from sklearn.mixture import GaussianMixture
from SimPEG import DataMisfit, Regularization, Optimization, InvProblem, Directives, Inversion
import SimPEG
import scipy.sparse as sp

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

def getCylinderPoints(xc,zc,r):
    xLocOrig1 = np.arange(-r,r+r/10.,r/10.)
    xLocOrig2 = np.arange(r,-r-r/10.,-r/10.)
    # Top half of cylinder
    zLoc1 = np.sqrt(-xLocOrig1**2.+r**2.)+zc
    # Bottom half of cylinder
    zLoc2 = -np.sqrt(-xLocOrig2**2.+r**2.)+zc
    # Shift from x = 0 to xc
    xLoc1 = xLocOrig1 + xc*np.ones_like(xLocOrig1)
    xLoc2 = xLocOrig2 + xc*np.ones_like(xLocOrig2)

    topHalf = np.vstack([xLoc1,zLoc1]).T
    topHalf = topHalf[0:-1,:]
    bottomHalf = np.vstack([xLoc2,zLoc2]).T
    bottomHalf = bottomHalf[0:-1,:]

    cylinderPoints = np.vstack([topHalf,bottomHalf])
    cylinderPoints = np.vstack([cylinderPoints,topHalf[0,:]])
    return cylinderPoints

cylinderPoints0 = getCylinderPoints(x0,z1,r0)
cylinderPoints1 = getCylinderPoints(x1,z1,r1)  

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
dmisAll = DataMisfit.l2_DataMisfit(survey)


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


from SimPEG.Maps import IdentityMap
from scipy.fftpack import dct,idct
class DCTMap(IdentityMap):
    """
        Changes the model into the physical property.

        If \\(p\\) is the physical property and \\(m\\) is the model, then

        .. math::

            p = \\log(m)

        and

        .. math::

            m = \\exp(p)

        NOTE: If you have a model which is log conductivity
        (ie. \\(m = \\log(\\sigma)\\)),
        you should be using an ExpMap

    """

    def __init__(self, mesh=None, nP=None, **kwargs):
        super(DCTMap, self).__init__(mesh=mesh, nP=nP, **kwargs)

    def _transform(self, m):
        return Utils.mkvc(dct(dct(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho'))

    def deriv(self, m, v=None):
        if v is not None:
            return dct(dct(v.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho')
        else:
            print "not implemented"

    def inverse(self, m):
        return Utils.mkvc(idct(idct(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho'))


class iDCTMap(IdentityMap):
    """
        Changes the physical proprety into the model

        If \\(p\\) is the physical property and \\(m\\) is the model, then

        .. math::

            p = \\log(m)

        and

        .. math::

            m = \\exp(p)

        NOTE: If you have a model which is log conductivity
        (ie. \\(m = \\log(\\sigma)\\)),
        you should be using an ExpMap

    """

    def __init__(self, mesh, nP=None, **kwargs):
        super(iDCTMap, self).__init__(mesh=mesh, nP=nP, **kwargs)

    def _transform(self, m):
        return Utils.mkvc(idct(idct(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho'))

    def deriv(self, m, v=None):
        if v is not None:
            return idct(idct(v.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho')
        else:
            print "not implemented"

    def inverse(self, m):
        return Utils.mkvc(dct(dct(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'), axis=0,norm = 'ortho'), axis=1,norm = 'ortho'))


import spgl1

#Parameter for SPGL1 iterations
nits = 10
mdct = (-5.)*np.ones_like(mtrue)
it = 0
#phi_d_normal = np.load('../phid_normal.npy')
#ratio = np.r_[6.5,phi_d_normal[0:-1]/phi_d_normal[1:]]
ratio = 10.*np.ones(nits)
min_progress = 1.2
xlist = []

#Parameters for W
nsubSrc = 5
InnerIt = 1
dmisfitsub = []
dmisfitall = []
dmisfitall.append(dmisAll.eval(mdct)/survey.nD)

#Initialize Random Source
W = np.random.randn(survey.nSrc,nsubSrc)
#problem.unpair()
#roblem.pair(survey)
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
dmis = DataMisfit.l2_DataMisfit(survey_r)
dmisfitsub.append(dmis.eval(mdct)/survey_r.nD)

problem.unpair()
problem.pair(survey)

print "end iteration: ",it, '; Overall Normalized Misfit: ', dmisAll.eval(mdct)/survey.nD

while (dmisAll.eval(mdct)/survey.nD)>0.5 and it<nits:
    
    problem.unpair()
    problem.pair(survey_r)

    def JS(x,mode):
        if mode == 1:
            return problem.Jvec(mdct,x)
        else:
            return problem.Jtvec(mdct,x)
    
    b = survey_r.dpred(mdct)-survey_r.dpred(mtrue)

    print "# of data: ", b.shape

    opts = spgl1.spgSetParms({'iterations':100, 'verbosity':2})
    sigtol = np.linalg.norm(b)/np.maximum(ratio[it],min_progress)
    #tautol = 20000.
    x,resid,grad,info = spgl1.spg_bpdn(JS, b, sigma = sigtol,options=opts)
    #x,resid,grad,info = spgl1.spg_lasso(JS,b,tautol,opts)
    #assert dmis.eval(mdct) > dmis.eval(mdct - x)
    mdct = mdct - x
    xlist.append(x)
    it +=1
    print "end iteration: ",it, '; Subsample Normalized Misfit: ', dmis.eval(mdct)/survey_r.nD
    dmisfitsub.append(dmis.eval(mdct)/survey_r.nD)

    problem.unpair()
    problem.pair(survey)
    dmisfitall.append(dmisAll.eval(mdct)/survey.nD)
    print "Dmisfit compared to full dataset: ",dmisAll.eval(mdct)/survey.nD

    if np.mod(it,InnerIt) ==0:
        W = np.random.randint(0, high=2, size=[survey.nSrc,nsubSrc])*2-1
        print 'update W'

        #problem.unpair()
        #roblem.pair(survey)
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
        dmis = DataMisfit.l2_DataMisfit(survey_r)

        d = survey_r.dpred(mtrue)
        survey_r.dobs = d
        survey_r.std = np.ones_like(d)*0.05
        survey_r.eps = 1e-5*np.linalg.norm(survey_r.dobs)
        print "end Update W; iteration: ",it, '; New Subsample Normalized Misfit: ', dmis.eval(mdct)/survey_r.nD

        problem.unpair()
        problem.pair(survey)
        
np.save('./dmisfitsub.npy',dmisfitsub)
np.save('./dmisfitall.npy',dmisfitall)
np.save('./mfinal.npy',mdct)
np.savez('./xlist.npz',xlist)

mm = mesh.plotImage(mdct)
plt.colorbar(mm[0])
plt.gca().set_xlim([-10.,10.])
plt.gca().set_ylim([-10.,0.])
plt.plot(cylinderPoints0[:,0],cylinderPoints0[:,1], linestyle = 'dashed', color='k')
plt.plot(cylinderPoints1[:,0],cylinderPoints1[:,1], linestyle = 'dashed', color='k')

plt.show()