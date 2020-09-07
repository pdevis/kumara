"""Run MCMC for atmospheric retrieval"""

'''___Built-In Modules___'''
from kumara.forwardModel.ForwardModel import ForwardModel
from kumara.forwardModel.UVspec import UVspec

'''___Third-Party Modules___'''
import os
import numpy as np
import matplotlib.pyplot as plt

'''___NPL Modules___'''


'''___Authorship___'''
__author__ = "Pieter De Vis"
__created__ = "11/11/2019"
__maintainer__ = "Pieter De Vis"
__email__ = "pieter.de.vis@npl.co.uk"
__status__ = "Development"



libradtranpath =  "/mnt/c/Users/annap/DATA/thesis/libRadtran-2.0.3/"

class LibradtranModel(ForwardModel):
    def __init__(self,path=os.getcwd()):
        """
        Initialise Forward Model

        :type product_path: str
        :param product_path: The data product file path.

        :type detail: str
        :param detail: Can take values:

        * "min" (default) - only information available with filename parsed.
        * "max" - information available in metadata files also parsed, but with opening the product.

        :type kwargs: -
        :param kwargs: Parsing parameters
        """

        # Initialise class attributes
        if not os.path.exists(path):
            os.makedirs(path)

        self.path=path
        self.lambd = None 
        self.TOAradiance = None
        self.RTcode = 'libradtran'

        self.edir = None  # Parse attributes dictionary
        self.edn = None  # Parse attributes dictionary
        self.eup = None  # Parse attributes dictionary
        self.uavg = None  # Parse attributes dictionary
        self.uu = None  # Parse attributes dictionary
        self.eglo = None  # Parse attributes dictionary

    def read_output(self,filename):
        self.lambd, self.edir, self.edn, self.eup, self.uavg, self.uu, self.eglo = np.genfromtxt('%s/%s.OUT'%(self.path,filename),unpack=True,usecols=[0,1,2,3,4,5,6])

  
    def run_model(self,filename,albedo_file,aerosol,H2O=-99,CH4=-99,O3=-99,CO2=-99,logvalues=False, verbose=False,sza=30,saa=180,vza=0,vaa=0,altitude=470,aerosol_type=1,atmosphere="midlatitude_summer",rte_solver='disort',wavs_range="400.0 2450.0",vert_squeeze=False):
        # some sensible alternative numbers for the variables:
        #H2Os = ["0","15","30"]
        #O3s = ["200","300","500"]
        #CH4s = [20e+18,30e+18,40e+18]
        #CO2s = [5e+21,7.5e+21,10e+21]
        #atmospheres = ["afglus","tropics","midlatitude_summer","midlatitude_winter","subarctic_summer","subarctic_winter","US-standard"]
        #solvers=['twostr','disort','MYSTIC']

        inputFilename = filename+'.INP'
        outputFilename = filename+'.OUT'
        inp = os.path.join(self.path,inputFilename)
        out = os.path.join(self.path,outputFilename)

        uvspec = UVspec()
        uvspec.inp["data_files_path"] = libradtranpath+'data'

        if verbose:
            uvspec.inp["verbose"]=""
        else:
            uvspec.inp["quiet"]=""    
        uvspec.inp["atmosphere_file"] = atmosphere
        uvspec.inp["source"] = 'solar '+libradtranpath+'data/solar_flux/kurudz_0.1nm.dat'


        if type(albedo_file)==str:
            uvspec.inp["albedo_file"] = albedo_file
        else:
            uvspec.inp["albedo_library IGBP"] = ""
            uvspec.inp["brdf_rpv_type"] = albedo_file

        #uvspec.inp["albedo_file"] = albedo_file
        #uvspec.inp["albedo"] = 0.5
        uvspec.inp["pseudospherical"] = ""

        uvspec.inp["sza"] = str(sza)
        uvspec.inp["umu"] = str(np.cos(np.deg2rad(vza)))
        uvspec.inp["phi"] = str(vaa)
        uvspec.inp["phi0"] = str(180-saa)

        uvspec.inp["altitude"] = str(altitude/1000)
        #if latitude>0:
            #uvspec.inp["latitude"] = "N "+str(latitude)
       # else:    
            #uvspec.inp["latitude"] = "S "+str(-latitude)
        #if longitude>0:
            #uvspec.inp["longitude"] = "E "+str(longitude)
        #else:    
            #uvspec.inp["longitude"] = "W "+str(-longitude)
        
        uvspec.inp["zout"] = "TOA"

        uvspec.inp["rte_solver"] = rte_solver

        uvspec.inp["wavelength"] = wavs_range  # Wavelength range [nm]
        uvspec.inp["mol_abs_param"] = 'reptran coarse'

        if H2O>-99:
            if logvalues==True:
                uvspec.inp["mol_modify H2O"] = str(10**H2O)+" MM"
            else:
                uvspec.inp["mol_modify H2O"] = str(H2O)+" MM"

        if CH4>-99:
            if logvalues==True:
                uvspec.inp["mol_modify CH4"] = str(10**CH4*10**18)+" cm_2"
            else:
                uvspec.inp["mol_modify CH4"] = str(CH4*10**18)+" cm_2"
        if O3>-99:
            if logvalues==True:
                uvspec.inp["mol_modify O3"] = str(10**O3)+" DU"
            else:
                uvspec.inp["mol_modify O3"] = str(O3)+" DU"
        if CO2>-99:
            if logvalues==True:
                uvspec.inp["mol_modify CO2"] = str(10**CO2)+" cm_2"
            else:
                uvspec.inp["mol_modify CO2"] = str(CO2)+" cm_2"

        #deltam off               # disable delta-scaling
        #uvspec.inp["number_of_streams"] = "24"      # number of streams used in DISORT

        uvspec.inp["aerosol_default"] = ""  # the simplest way to include aerosol :-)
        #uvspec.inp["aerosol_modify tau scale"] = aerosol 
        uvspec.inp["aerosol_set_tau_at_wvl 440"] = aerosol


        # uvspec.inp["aerosol_vulcan"] = 1          # Aerosol type above 2km
        uvspec.inp["aerosol_haze"] =  aerosol_type            # Aerosol type below 2km
        #uvspec.inp["aerosol_haze"] =  1            # Aerosol type below 2km
        # uvspec.inp["aerosol_season"] =  1          # Summer season
        # uvspec.inp["aerosol_visibility"] =  50.0   # 
        
        # #uvspec.inp["aerosol_angstrom"] = "1.1 0.07" # Scale aerosol optical depth 
        #                   # using Angstrom alpha and beta
        #                   # coefficients
        # #uvspec.inp["aerosol_modify ssa scale"] =  0.85    # Scale the single scattering albedo 
        #                   # for all wavelengths
        # #uvspec.inp["aerosol_modify gg set"] =  0.70       # Set the asymmetry factor
        # uvspec.inp["aerosol_file tau"] =  os.path.join(self.path,"AERO_TAU.DAT")
        # uvspec.inp["aerosol_sizedist_file"] =  os.path.join(self.path,"AERO_SIZEDIST.DAT")
        
        # uvspec.inp["aerosol_refrac_index "] = "1.55 0.002"
        # #uvspec.inp["aerosol_refrac_file "] = os.path.join(self.path,"AERO_REFRAC.DAT")
        #                   # File with aerosol refractive index
        #uvspec.inp["disort_intcor"] = "moments"

        #uvspec.inp["profile_file aer1 1D"] = os.path.join(self.path,"testwc.DAT") + " "
        #uvspec.inp["profile_properties aer1"] = os.path.join(self.path,"aerosolmie.cdf") + " interpolate"
        #if logvalues:
            #uvspec.inp["profile_modify aer1 tau set"] = str(10**aerosol)
        #else:
            #uvspec.inp["profile_modify aer1 tau set"] = str(aerosol)

        if vert_squeeze:
            uvspec.inp['aerosol_profile_modtran']=''


        uvspec.inp["output_user"] = "lambda edir edn eup uavg uu eglo"
        #uvspec.inp["output_quantity"] = 'reflectivity' #'transmittance' #

        uvspec.write_input(inp)
        uvspec.run(inp,out,verbose,path=libradtranpath)
        return None

    def get_TOA(self,filename):
        rad = np.genfromtxt('%s/%s.OUT'%(self.path,filename),unpack=True,usecols=[5]) #in mW / (m2 nm)
        rad[np.isnan(rad)]=0
        return rad

    def get_wavs(self,filename):
        wavs = np.genfromtxt('%s/%s.OUT'%(self.path,filename),unpack=True,usecols=[0])
        return wavs     

    def plot_spectrum(self,filename):
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        ax.plot(self.lambd,self.edir)
        ax.plot(self.lambd,self.edn)
        ax.plot(self.lambd,self.eup)
        ax.plot(self.lambd,self.uavg)
        ax.plot(self.lambd,self.uu)
        ax.plot(self.lambd,self.eglo)
        ax.set_xlabel(r"$\lambda$ (nm)")
        ax.set_ylabel(r"Brightness")
        ax.set_ylim([0,1])
        fig.savefig('%s/%s.png'%(self.path,filename))
        del fig
        return None
