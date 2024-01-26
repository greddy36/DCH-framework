import CRABClient
from CRABClient.UserUtilities import config
from WMCore.Configuration import Configuration

wdir='/uscms_data/d3/alkaloge/MetStudies/nAOD/CMSSW_10_6_5/src/newCrab/'

config = config()
config.General.workArea = 'Crab_projects_JECs'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.psetName = '{0:s}/Files/pset.py'.format(wdir)
#config.JobType.pluginName = 'PrivateMC'
config.JobType.inputFiles =[ '{0:s}/Files/make_jmev2.py'.format(wdir), '{0:s}/Files/keep_and_drop.txt'.format(wdir),  '{0:s}/Files/pset.py'.format(wdir), '{0:s}/Files/FrameworkJobReport.xml'.format(wdir)]
#config.JobType.outputFiles = ['all_HWminus_002_4of10.root']
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/lpcsusyhiggs/ntuples/nAODv9/'
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
DATASETHERE

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    #from WMCore.Configuration import Configuration
    #from CRABClient.UserUtilities import config

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_multiPub'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
