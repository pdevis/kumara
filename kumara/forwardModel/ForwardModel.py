"""Run MCMC for atmospheric retrieval"""

'''___Built-In Modules___'''

'''___Third-Party Modules___'''
#import Libradtran
import os
from abc import ABC, abstractmethod
'''___NPL Modules___'''


'''___Authorship___'''
__author__ = "Pieter De Vis"
__created__ = "11/11/2019"
__maintainer__ = "Pieter De Vis"
__email__ = "pieter.de.vis@npl.co.uk"
__status__ = "Development"

class ForwardModel(ABC): 
    def __init__(self, RTcode='libradtran',path='./libradtran/'):
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
        self.RTcode = RTcode

    

    @abstractmethod
    def read_output(self,filename):
        pass

    @abstractmethod
    def get_TOA(self,filename):
        pass

    @abstractmethod
    def get_wavs(self,filename):
        pass

    @abstractmethod
    def run_model(self,filename,albedo_file,aerosol,**kwargs):
        pass

    @abstractmethod
    def plot_spectrum(self,filename):
        pass   
