#!/usr/bin/python3

import argparse
import sys
import json
import pprint
import re
import time
import traceback
import shlex

def error(msg):
    sys.stderr.write(msg)
    sys.exit()

def asp_to_loesung(path) -> str:
    f = open(path, 'r')

    content = f.readlines()
    if len(content) == 0:
        error('error - ' + path + ' is empty\n')

    answer = ""
    for i in content:
        if 'error' in i or 'ERROR' in i:
            return error("error found in solution")
        if "dl" in i:
            answer = i
    answer = answer.split(' ')
    answer = [item for item in answer if len(item) > 1]
    out = {}
    out["hash"] = 1 # no way to precompute this
    out["zugfahrten"] = []
    reghash    = re.compile('hash\((-?\d+)\)')
    regroute   = re.compile('route\("(.*)",(\d+),(.*)\)')
    regstart   = re.compile('dl\(\("(.*)",(\(.*\))\),"(\d+)"\)$')
    regfahrweg = re.compile('fahrweg\("(.*)","(.*)"\)')
    regtop     = re.compile('topologischeReihenfolge\("(.*)",(\d+),(\d+)\)')
    regkenn    = re.compile('abschnittskennzeichen\("(.*)",(\d+),"(.*)"\)')
    regfolge   = re.compile('abschnitt\(folgenid,"(.*)",(\d+),"(.*)"\)')
    redge      = re.compile('edge\("(.*)",(\d+),(\(.*\)),(\(.*\))\)$')
    visit    = {}
    times    = {}
    result   = {}
    sifw     = {}
    topology = {}
    kenn     = {}
    folge    = {}
    ingoing  = {}
    outgoing = {}

    for i in answer:
        x = regstart.match(i)
        if x != None:
            times.setdefault(x.group(1),{})[str(x.group(2))] = int(x.group(3)) 
        x = reghash.match(i)
        if x != None:
            out["verkehrsplanHash"] = x.group(1)
        x = regroute.match(i)
        if x != None:
            visit.setdefault(x.group(1),[]).extend([x.group(2),x.group(3)])
        x = regfahrweg.match(i)
        if x != None:
            sifw[x.group(1)] = x.group(2)
        x = regtop.match(i)
        if x != None:
            if x.group(1) not in topology:
                topology[x.group(1)] = {}
            topology[x.group(1)][x.group(2)] = int(x.group(3))
        x = regkenn.match(i)
        if x != None:
            if x.group(1) not in kenn:
                kenn[x.group(1)] = {}
            kenn[x.group(1)][x.group(2)] = x.group(3)
        x = regfolge.match(i)
        if x != None:
            if x.group(1) not in folge:
                folge[x.group(1)] = {}
            folge[x.group(1)][x.group(2)] = x.group(3)
        x = redge.match(i)
        if x != None:
            ingoing.setdefault(x.group(1),{}).setdefault(x.group(3),[]).append(x.group(2))
            outgoing.setdefault(x.group(1),{}).setdefault(x.group(4),[]).append(x.group(2))

    for si in times:
        for v,t in times[si].items():
            if v in ingoing[si]:
                for i in ingoing[si][v]:
                    if i in visit[si]:
                        result.setdefault(si,{}).setdefault(i,{})["ein"] = str(time.strftime('%H:%M:%S', time.gmtime(t)))
            if v in outgoing[si]:
                for i in outgoing[si][v]:
                    if i in visit[si]:
                        result.setdefault(si,{}).setdefault(i,{})["aus"] = str(time.strftime('%H:%M:%S', time.gmtime(t))) 

    ins = out["zugfahrten"]
    for i in result:
        si = {}
        si["funktionaleAngebotsbeschreibungId"]=i
        si["zugfahrtabschnitte"]=[]
        zfa = si["zugfahrtabschnitte"]
        for k in result[i]:
            abschnitt={}
            abschnitt["fahrwegabschnittId"]=sifw[i]+"#"+str(str(k))
            abschnitt["fahrweg"]=sifw[i]
            abschnitt["ein"]=result[i][k]["ein"]
            abschnitt["aus"]=result[i][k]["aus"]
            abschnitt["abschnittsfolge"]=folge[i][k]
            abschnitt["reihenfolge"] = topology[i][k]
            if i in kenn and k in kenn[i]:
                abschnitt["abschnittsvorgabe"] = kenn[i][k]
            zfa.append(abschnitt)

        ins.append(si)
    json.dump(out, sys.stdout, indent=4, sort_keys=True)
    f.close()
    return None

def main():
    try:
        parser = argparse.ArgumentParser(description="Converts a clingoDL answer set into a challenge solution ")
        parser.add_argument('input', nargs='?', metavar='<input>', type=argparse.FileType('r'), default=sys.stdin, help='Read from %(metavar)s (default: <stdin>)')

        args = parser.parse_args()
        asp_to_loesung(args.input.name)
    except Exception as e:
        traceback.print_exception(*sys.exc_info())
        return 1

if __name__ == '__main__':
    sys.exit(main())
