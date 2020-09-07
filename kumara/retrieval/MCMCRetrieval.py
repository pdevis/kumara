"""Run MCMC for atmospheric retrieval"""

'''___Built-In Modules___'''
from kumara.forwardModel.ForwardModelFactory import ForwardModelFactory
from kumara.sensors.SensorFactory import SensorFactory

'''___Third-Party Modules___'''
import numpy as np
import emcee
import threading
from multiprocessing import Pool
import time
import os
from scipy.linalg import lapack

os.environ["OMP_NUM_THREADS"] = "1"


'''___NPL Modules___'''


'''___Authorship___'''
__author__ = "Pieter De Vis"
__created__ = "01/03/2020"
__maintainer__ = "Pieter De Vis"
__email__ = "pieter.de.vis@npl.co.uk"
__status__ = "Development"


inds_cache = {}


# lock = threading.Lock()
class MCMCRetrieval:
    def __init__(self,sensor,filename,wavs,FWHM,observed,reflpath,altitude,wavs_range,aerosol_type,vert_squeeze,sza,vza,vaa,saa,uncertainty=None,cov=None,logvalues=False,uplims=+np.inf,downlims=-np.inf,forwardmodel='libradtran',path='MCMC'):
        self.filename=filename
        self.observed=observed
        self.wavs=wavs
        self.FWHM=FWHM
        self.uncertainty=uncertainty
        self.cov=cov
        self.FM=ForwardModelFactory.create_model(forwardmodel)
        self.sensor=SensorFactory.create_sensor(sensor,wavs,FWHM)
        self.reflpath=reflpath
        self.uplims=np.array(uplims)
        self.downlims=np.array(downlims)
        self.logvalues=logvalues
        self.altitude=altitude
        #self.latitude=latitude
        #self.longitude=longitude
        self.sza=sza
        self.vza=vza
        self.vaa=vaa
        self.saa=saa
        self.aerosol_type=aerosol_type
        self.wavs_range=wavs_range
        self.vert_squeeze=vert_squeeze

    def run_retrieval(self,theta_0,nwalkers,steps):
        
        ndimw = len(theta_0)
        pos = [theta_0*np.random.normal(1.0,0.1,len(theta_0)) for i in range(nwalkers)]
        # print(self.lnprob(theta_0))
        with Pool(processes=20) as pool:
           sampler = emcee.EnsembleSampler(nwalkers, ndimw, self.lnprob,pool=pool)
           sampler.run_mcmc(pos, steps, progress=True)
        # sampler = emcee.EnsembleSampler(nwalkers, ndimw, self.lnprob)
        # sampler.run_mcmc(pos, steps, progress=True)
        return sampler.get_chain()[:].reshape(-1,ndimw)

    def find_chisum(self,theta):
        # lock.acquire()
        filenameMCMC=self.filename+"_MCMC_%s"%time.time()
        # lock.release()

        self.FM.run_model(filenameMCMC,self.reflpath,*theta,logvalues=self.logvalues,altitude=self.altitude, aerosol_type=self.aerosol_type,wavs_range=self.wavs_range,vert_squeeze=self.vert_squeeze,sza=self.sza, vza=self.vza, vaa=self.vaa, saa=self.saa)
        rad_RT=self.FM.get_TOA(filenameMCMC)
        wavs_RT=self.FM.get_wavs(filenameMCMC)
        model=self.sensor.convolve(wavs_RT,rad_RT)
        diff=model-self.observed
        if self.cov is None:
            return np.sum((diff)**2/self.uncertainty**2)
        else:
            cov=np.ascontiguousarray(self.cov)
            return np.abs(np.dot(np.dot(diff.T,np.linalg.inv(cov)),diff))

    # def upper_triangular_to_symmetric(self,ut):
    #     n = ut.shape[0]
    #     try:
    #         inds = inds_cache[n]
    #     except KeyError:
    #         inds = np.tri(n, k=-1, dtype=np.bool)
    #         inds_cache[n] = inds
    #     ut[inds] = ut.T[inds]


    # def fast_positive_definite_inverse(self,m):
    #     cholesky, info = lapack.dpotrf(m)
    #     if info != 0:
    #         raise ValueError('dpotrf failed on input {}'.format(m))
    #     inv, info = lapack.dpotri(cholesky)
    #     if info != 0:
    #         raise ValueError('dpotri failed on input {}'.format(cholesky))
    #     self.upper_triangular_to_symmetric(inv)
    #     return inv


    def lnlike(self,theta):
        #print(theta,[10**theta[0],10**theta[1]],self.find_chisum(theta))
        return -0.5*(self.find_chisum(theta))

    def lnprior(self,theta):
        if all(self.downlims<theta) and all(self.uplims>theta):
            return 0.0
        else:    
            return -np.inf

    def lnprob(self,theta):
        lp_prior = self.lnprior(theta)
        if not np.isfinite(lp_prior):
            return -np.inf
        lp=self.lnlike(theta)
        #if lp<-999:
            #return -np.inf
        if not np.isfinite(lp):
            return -np.inf    
        return lp_prior + lp
