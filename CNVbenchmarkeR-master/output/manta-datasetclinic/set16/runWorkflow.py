#!/usr/bin/env python2
# Workflow run script auto-generated by command: '/home/robin/Documents/Project/manta/Install/bin/configManta.py --bam=/home/robin/Documents/Project/Samples/bam/all/17398.bam --bam=/home/robin/Documents/Project/Samples/bam/all/17399.bam --bam=/home/robin/Documents/Project/Samples/bam/all/17400.bam --bam=/home/robin/Documents/Project/Samples/bam/all/17401.bam --bam=/home/robin/Documents/Project/Samples/bam/all/17402.bam --bam=/home/robin/Documents/Project/Samples/bam/all/17403.bam --referenceFasta=/home/robin/Documents/Project/Samples/hg38/hg38.fa --config=/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/manta-datasetclinic/configManta.py.ini --exome --runDir=/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/manta-datasetclinic/set16'
#

import os, sys

if sys.version_info >= (3,0):
    import platform
    raise Exception("Manta does not currently support python3 (version %s detected)" % (platform.python_version()))

if sys.version_info < (2,6):
    import platform
    raise Exception("Manta requires python2 version 2.6+ (version %s detected)" % (platform.python_version()))

scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(r'/home/robin/Documents/Project/manta/Install/lib/python')

from mantaWorkflow import MantaWorkflow



def get_run_options(workflowClassName) :

    from optparse import OptionGroup, SUPPRESS_HELP

    from configBuildTimeInfo import workflowVersion
    from configureUtil import EpilogOptionParser
    from estimateHardware import EstException, getNodeHyperthreadCoreCount, getNodeMemMb

    epilog="""Note this script can be re-run to continue the workflow run in case of interruption.
Also note that dryRun option has limited utility when task definition depends on upstream task
results -- in this case the dry run will not cover the full 'live' run task set."""

    parser = EpilogOptionParser(description="Version: %s" % (workflowVersion), epilog=epilog, version=workflowVersion)


    parser.add_option("-m", "--mode", type="string",dest="mode",
                      help=SUPPRESS_HELP)
    parser.add_option("-j", "--jobs", type="string",dest="jobs",
                  help="number of jobs, must be an integer or 'unlimited' (default: Estimate total cores on this node)")
    parser.add_option("-g","--memGb", type="string",dest="memGb",
                  help="gigabytes of memory available to run workflow, must be an integer (default: Estimate the total memory for this node)")
    parser.add_option("-d","--dryRun", dest="isDryRun",action="store_true",default=False,
                      help="dryRun workflow code without actually running command-tasks")
    parser.add_option("--quiet", dest="isQuiet",action="store_true",default=False,
                      help="Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)")

    def isLocalSmtp() :
        import smtplib
        try :
            smtplib.SMTP('localhost')
        except :
            return False
        return True

    isEmail = isLocalSmtp()
    emailHelp = SUPPRESS_HELP
    if isEmail :
        emailHelp="send email notification of job completion status to this address (may be provided multiple times for more than one email address)"

    parser.add_option("-e","--mailTo", type="string",dest="mailTo",action="append",help=emailHelp)

    debug_group = OptionGroup(parser,"development debug options")
    debug_group.add_option("--rescore", dest="isRescore",action="store_true",default=False,
                          help="Reset task list to re-run hypothesis generation and scoring without resetting graph generation.")

    parser.add_option_group(debug_group)

    ext_group = OptionGroup(parser,"extended portability options (should not be needed by most users)")
    ext_group.add_option("--maxTaskRuntime", type="string", metavar="hh:mm:ss",
                      help="Specify max runtime per task (no default)")

    parser.add_option_group(ext_group)

    (options,args) = parser.parse_args()

    if not isEmail : options.mailTo = None

    if len(args) :
        parser.print_help()
        sys.exit(2)

    if options.mode is None :
        options.mode = "local"
    elif options.mode not in ["local"] :
        parser.error("Invalid mode. Available modes are: local")

    if options.jobs is None :
        try :
            options.jobs = getNodeHyperthreadCoreCount()
        except EstException:
            parser.error("Failed to estimate cores on this node. Please provide job count argument (-j).")
    if options.jobs != "unlimited" :
        options.jobs=int(options.jobs)
        if options.jobs <= 0 :
            parser.error("Jobs must be 'unlimited' or an integer greater than 1")

    # note that the user sees gigs, but we set megs
    if options.memGb is None :
        try :
            options.memMb = getNodeMemMb()
        except EstException:
            parser.error("Failed to estimate available memory on this node. Please provide available gigabyte argument (-g).")
    elif options.memGb != "unlimited" :
        options.memGb=int(options.memGb)
        if options.memGb <= 0 :
            parser.error("memGb must be 'unlimited' or an integer greater than 1")
        options.memMb = 1024*options.memGb
    else :
        options.memMb = options.memGb

    options.resetTasks=[]
    if options.isRescore :
        options.resetTasks.append("makeHyGenDir")

    return options


def main(pickleConfigFile, primaryConfigSection, workflowClassName) :

    from configureUtil import getConfigWithPrimaryOptions

    runOptions=get_run_options(workflowClassName)
    flowOptions,configSections=getConfigWithPrimaryOptions(pickleConfigFile,primaryConfigSection)

    # new logs and marker files to assist automated workflow monitoring:
    warningpath=os.path.join(flowOptions.runDir,"workflow.warning.log.txt")
    errorpath=os.path.join(flowOptions.runDir,"workflow.error.log.txt")
    exitpath=os.path.join(flowOptions.runDir,"workflow.exitcode.txt")

    # the exit path should only exist once the workflow completes:
    if os.path.exists(exitpath) :
        if not os.path.isfile(exitpath) :
            raise Exception("Unexpected filesystem item: '%s'" % (exitpath))
        os.unlink(exitpath)

    wflow = workflowClassName(flowOptions)

    retval=1
    try:
        retval=wflow.run(mode=runOptions.mode,
                         nCores=runOptions.jobs,
                         memMb=runOptions.memMb,
                         dataDirRoot=flowOptions.workDir,
                         mailTo=runOptions.mailTo,
                         isContinue="Auto",
                         isForceContinue=True,
                         isDryRun=runOptions.isDryRun,
                         isQuiet=runOptions.isQuiet,
                         resetTasks=runOptions.resetTasks,
                         successMsg=wflow.getSuccessMessage(),
                         retryWindow=0,
                         retryMode='all',
                         warningLogFile=warningpath,
                         errorLogFile=errorpath)
    finally:
        exitfp=open(exitpath,"w")
        exitfp.write("%i\n" % (retval))
        exitfp.close()

    sys.exit(retval)

main(r"/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/manta-datasetclinic/set16/runWorkflow.py.config.pickle","manta",MantaWorkflow)

