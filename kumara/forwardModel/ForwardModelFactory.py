"""Run MCMC for atmospheric retrieval"""

'''___Built-In Modules___'''
from kumara.forwardModel.LibradtranModel import LibradtranModel

'''___Third-Party Modules___'''


'''___NPL Modules___'''


'''___Authorship___'''
__author__ = "Pieter De Vis"
__created__ = "11/11/2019"
__maintainer__ = "Pieter De Vis"
__email__ = "pieter.de.vis@npl.co.uk"
__status__ = "Development"


class ForwardModelFactory():
    @staticmethod
    def create_model(RTcode,*args,**kwargs):
        if RTcode=='libradtran':
            print("Using libradtran as forward model. \n")
            return LibradtranModel(*args,**kwargs)


