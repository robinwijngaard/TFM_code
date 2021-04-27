#!/usr/bin/env python
#
# This is a script to create a dot graph from pyflow state files.
# Usage: $script >| task_graph.dot
#
# Note that script assumes the default pyflow state files are in the script directory.
#
# This file was autogenerated by process: '/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master'
# ...from working directory: '/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/manta-datasetall/set16/runWorkflow.py'
#

import datetime,os,sys,time

scriptDir=os.path.abspath(os.path.dirname(__file__))


def timeStampToTimeStr(ts) :
    """
    converts time.time() output to timenow() string
    """
    return datetime.datetime.utcfromtimestamp(ts).isoformat()+"Z"


def timeStrNow():
    return timeStampToTimeStr(time.time())


def cmdline() :
    return " ".join(sys.argv)


class Bunch:
    """
    generic struct with named argument constructor
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class LogGlobals :
    isFsync = True


def hardFlush(ofp):
    ofp.flush()
    if ofp.isatty() : return
    # fsync call has been reported to consistently fail in some contexts (rsh?)
    # so allow OSError
    if not LogGlobals.isFsync : return
    try :
        os.fsync(ofp.fileno())
    except OSError:
        LogGlobals.isFsync = False


class TaskNodeConstants(object) :

    validRunstates = ("complete", "running", "queued", "waiting", "error")


class DotConfig(object) :
    """
    A static container of configuration data for dot graph output
    """

    runstateDotColor = {"waiting" : "grey",
                        "running" : "green",
                        "queued" : "yellow",
                        "error" : "red",
                        "complete" : "blue" }

    runstateDotStyle = {"waiting" : "dashed",
                        "running" : None,
                        "queued" : None,
                        "error" : "bold",
                        "complete" : None }

    @staticmethod
    def getRunstateDotAttrib(runstate) :
        color = DotConfig.runstateDotColor[runstate]
        style = DotConfig.runstateDotStyle[runstate]
        attrib = ""
        if color is not None : attrib += " color=%s" % (color)
        if style is not None : attrib += " style=%s" % (style)
        return attrib

    @staticmethod
    def getTypeDotAttrib(nodeType) :
        attrib = ""
        if nodeType == "workflow" :
            attrib += " shape=rect style=rounded"
        return attrib

    @staticmethod
    def getDotLegend() :
        string = '{ rank = source; Legend [shape=none, margin=0, label=<\n'
        string += '<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">\n'
        string += '<TR><TD COLSPAN="2">Legend</TD></TR>\n'
        for state in TaskNodeConstants.validRunstates :
            color = DotConfig.runstateDotColor[state]
            string += '<TR> <TD>%s</TD> <TD BGCOLOR="%s"></TD> </TR>\n' % (state, color)
        string += '</TABLE>>];}\n'
        return string


def taskStateParser(stateFile) :
    class Constants :
        nStateCols = 5

    for line in open(stateFile) :
        if len(line) and line[0] == "#" : continue
        line = line.strip()
        w = line.split("\t")
        if len(w) != Constants.nStateCols :
            raise Exception("Unexpected format in taskStateFile: '%s' line: '%s'" % (stateFile, line))
        yield [x.strip() for x in w]


def taskInfoParser(infoFile) :
    class Constants :
        nInfoCols = 10

    for line in open(infoFile) :
        if len(line) and line[0] == "#" : continue
        line = line.lstrip()
        w = line.split("\t", (Constants.nInfoCols - 1))
        if len(w) != Constants.nInfoCols :
            raise Exception("Unexpected format in taskInfoFile: '%s' line: '%s'" % (infoFile, line))
        yield [x.strip() for x in w]


def getTaskInfoDepSet(s) :
    # reconstruct dependencies allowing for extraneous whitespace in the file:
    s = s.strip()
    if s == "" : return []
    return set([d.strip() for d in s.split(",")])


def writeDotGraph(taskInfoFile, taskStateFile, workflowClassName) :
    """
    write out the current graph state in dot format
    """

    addOrder = []
    taskInfo = {}
    headNodes = set()
    tailNodes = set()

    # read info file:
    for (label, namespace, ptype, _nCores, _memMb, _priority, _isForceLocal, depStr, _cwdStr, _command) in taskInfoParser(taskInfoFile) :
        tid = (namespace, label)
        addOrder.append(tid)
        taskInfo[tid] = Bunch(ptype=ptype,
                              parentLabels=getTaskInfoDepSet(depStr))
        if len(taskInfo[tid].parentLabels) == 0 : headNodes.add(tid)
        tailNodes.add(tid)
        for plabel in taskInfo[tid].parentLabels :
            ptid = (namespace, plabel)
            if ptid in tailNodes : tailNodes.remove(ptid)

    for (label, namespace, runState, _errorCode, _time) in taskStateParser(taskStateFile) :
        tid = (namespace, label)
        taskInfo[tid].runState = runState

    dotFp = sys.stdout
    dotFp.write("// Task graph from pyflow object '%s'\n" % (workflowClassName))
    dotFp.write("// Process command: '%s'\n" % (cmdline()))
    dotFp.write("// Process working dir: '%s'\n" % (os.getcwd()))
    dotFp.write("// Graph capture time: %s\n" % (timeStrNow()))
    dotFp.write("\n")
    dotFp.write("digraph %s {\n" % (workflowClassName + "Graph"))
    dotFp.write("\tcompound=true;\nrankdir=LR;\nnode[fontsize=10];\n")
    labelToSym = {}
    namespaceGraph = {}
    for (i, (namespace, label)) in enumerate(addOrder) :
        tid = (namespace, label)
        if namespace not in namespaceGraph :
            namespaceGraph[namespace] = ""
        sym = "n%i" % i
        labelToSym[tid] = sym
        attrib1 = DotConfig.getRunstateDotAttrib(taskInfo[tid].runState)
        attrib2 = DotConfig.getTypeDotAttrib(taskInfo[tid].ptype)
        namespaceGraph[namespace] += "\t\t%s [label=\"%s\"%s%s];\n" % (sym, label, attrib1, attrib2)

    for (namespace, label) in addOrder :
        tid = (namespace, label)
        sym = labelToSym[tid]
        for plabel in taskInfo[tid].parentLabels :
            ptid = (namespace, plabel)
            namespaceGraph[namespace] += ("\t\t%s -> %s;\n" % (labelToSym[ptid], sym))

    for (i, ns) in enumerate(namespaceGraph.keys()) :
        isNs = ((ns is not None) and (ns != ""))
        dotFp.write("\tsubgraph cluster_sg%i {\n" % (i))
        if isNs :
            dotFp.write("\t\tlabel = \"%s\";\n" % (ns))
        else :
            dotFp.write("\t\tlabel = \"%s\";\n" % (workflowClassName))
        dotFp.write(namespaceGraph[ns])
        dotFp.write("\t\tbegin%i [label=\"begin\" shape=diamond];\n" % (i))
        dotFp.write("\t\tend%i [label=\"end\" shape=diamond];\n" % (i))
        for (namespace, label) in headNodes :
            if namespace != ns : continue
            sym = labelToSym[(namespace, label)]
            dotFp.write("\t\tbegin%i -> %s;\n" % (i, sym))
        for (namespace, label) in tailNodes :
            if namespace != ns : continue
            sym = labelToSym[(namespace, label)]
            dotFp.write("\t\t%s -> end%i;\n" % (sym, i))
        dotFp.write("\t}\n")
        if ns in labelToSym :
            dotFp.write("\t%s -> begin%i [style=dotted];\n" % (labelToSym[ns], i))
            # in LR orientation this will make the graph look messy:
            # dotFp.write("\tend%i -> %s [style=invis];\n" % (i,labelToSym[ns]))

    dotFp.write(DotConfig.getDotLegend())
    dotFp.write("}\n")
    hardFlush(dotFp)


if __name__ == '__main__' :
    writeDotGraph(os.path.join(scriptDir,'pyflow_tasks_info.txt'),os.path.join(scriptDir,'pyflow_tasks_runstate.txt'),'MantaWorkflow')

