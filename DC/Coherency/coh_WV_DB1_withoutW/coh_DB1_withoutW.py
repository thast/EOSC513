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

#2D model
csx, csy, csz = 0.25,0.25,0.25
# Number of core cells in each directiPon s
ncx, ncz = 2**7-24,2**7-12
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
dmis = DataMisfit.l2_DataMisfit(survey)
survey.dpred(mtrue)
survey.makeSyntheticData(mtrue,std=0.05,force=True)
survey.eps = 1e-5*np.linalg.norm(survey.dobs)

print '# of data: ', survey.dobs.shape

from SimPEG.Maps import IdentityMap
import pywt

class WaveletMap(IdentityMap):

    def __init__(self, mesh=None, nP=None, **kwargs):
        super(WaveletMap, self).__init__(mesh=mesh, nP=nP, **kwargs)

    def _transform(self, m, wv = 'db1'):
        coeff_wv = pywt.wavedecn(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
        array_wv = pywt.coeffs_to_array(coeff_wv)        
        return Utils.mkvc(array_wv[0])

    def deriv(self, m, v=None, wv = 'db1'):
        if v is not None:
            coeff_wv = pywt.wavedecn(v.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
            array_wv = pywt.coeffs_to_array(coeff_wv)        
            return Utils.mkvc(array_wv[0])
        else:
            print "not implemented"

    def inverse(self, m, wv = 'db1'):
        msyn = np.zeros(mesh.nC)
        coeff_wv = pywt.wavedecn(msyn.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
        array_wv = pywt.coeffs_to_array(coeff_wv)
        coeff_back = pywt.array_to_coeffs(m.reshape(array_wv[0].shape, order = 'F'),array_wv[1])
        coeff_m = pywt.waverecn(coeff_back,wv, mode = 'per')
        return Utils.mkvc(coeff_m)

class iWaveletMap(IdentityMap):

    def __init__(self, mesh, nP=None, **kwargs):
        super(iWaveletMap, self).__init__(mesh=mesh, nP=nP, **kwargs)

    def _transform(self, m, wv = 'db1'):
        msyn = np.zeros(mesh.nC)
        coeff_map = pywt.wavedecn(msyn.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
        array_map = pywt.coeffs_to_array(coeff_map)
        coeff_map = pywt.array_to_coeffs(m.reshape(array_map[0].shape,order= 'F'),array_map[1])
        coeff_back_map = pywt.waverecn(coeff_map,wv, mode = 'per')
        return Utils.mkvc(coeff_back_map)

    def deriv(self, m, v=None, wv = 'db1'):
        if v is not None:
            coeff_wv = pywt.wavedecn(v.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
            array_wv = pywt.coeffs_to_array(coeff_wv)
            coeff_back = pywt.array_to_coeffs(v,array_wv[1])
            coeff_m = pywt.waverecn(coeff_back,wv, mode = 'per')
            return Utils.mkvc(coeff_m)        
        else:
            print "not implemented"

    def inverse(self, m, wv = 'db1'):
        
        coeff_wv = pywt.wavedecn(m.reshape(self.mesh.nCx,self.mesh.nCy,order = 'F'),wv, mode = 'per')
        array_wv = pywt.coeffs_to_array(coeff_wv)        
        return Utils.mkvc(array_wv[0])

wavmap = WaveletMap(mesh)
iwavmap = iWaveletMap(mesh) 

J = lambda v:  problem.Jvec(mtrue,v)
Jt = lambda v: problem.Jtvec(mtrue,v)

x = np.zeros_like(mtrue)
v = np.zeros(problem.survey.nD)
print 'v shape: ',v.shape
indx = np.random.permutation(len(x))
indv = np.random.permutation(len(v))

coeff = 50
muwav = np.zeros([coeff, coeff])
for i in range(coeff):
    print 'iteration: ', i
    coeffs = pywt.wavedec2(np.zeros((mesh.nCx,mesh.nCy)),'db1')
    indx0 = np.random.randint(1,len(coeffs))
    indx1 = np.random.randint(0,3)
    indx2 = np.random.randint(0,(coeffs[indx0][indx1]).shape[0])
    indx3 = np.random.randint(0,(coeffs[indx0][indx1]).shape[1])
    coeffs[indx0][indx1][indx2][indx3] = 1.
    for j in range(coeff):
        v = np.zeros(problem.survey.nD)
        v[indv[j]] = 1.
        v = Jt(v)
        v = v/np.linalg.norm(v)
        x = pywt.waverec2(coeffs,'db1')
        x = x.reshape(mesh.nC,order = 'F')
        x = x/np.linalg.norm(x)
        
        muwav [i,j] = x.dot(v)


np.save('./mu.npy',muwav)
np.savez('./mu.npz',muwav)

fig1 = plt.figure(figsize=(10,8))
ax1 = fig1.add_subplot(111)
mm = ax1.matshow(np.log10(np.abs(muwav)), cmap="jet")
mm.set_clim(vmin=-8., vmax=0.)
cb = plt.colorbar(mm)
ax1.set_aspect("equal")
cb.set_label("Log10 |< * , * >|", fontsize=20)
cb.ax.tick_params(labelsize=20)
ax1.set_xlabel("columns of $J^T$",fontsize = 20)
ax1.set_ylabel("columns of $S^H$",fontsize =20)
ax1.set_title("Coherency of $J^T$ with $S^H$, $S$ the Wavelets Transform DB1",fontsize = 22)
ax1.text(0.,55.,"Minimum = %g"%(np.abs(muwav).min()),fontsize =20)
ax1.text(0.,57.,"Maximum = %.2f"%(np.abs(muwav).max()),fontsize =20)
ax1.tick_params(labelsize = 16)

fig1.savefig('./coh_DB1_withoutW.png')

plt.show()