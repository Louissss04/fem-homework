#!/usr/bin/python
# 
# This is the preprocessor (reader) for the script 'quadrule-from-csv'.
# It tries to identify and convert symmetry orbits in the input data.
#
# $Id: quadrule-from-csv.py,v 1.13 2020/03/18 10:12:52 zlb Exp $

from __future__ import print_function
import sys
import os
import re

###############################################################################
# define functions

# Use mpmath (pip install mpmath) for high-precision arithmetics
import mpmath
from mpmath import mp

if os.getenv('QUADRULE_DIGITS') == None:
    mp.dps = 34
else:
    mp.dps = int(os.environ['QUADRULE_DIGITS'])

# for equality test
if os.getenv('QUADRULE_EPS') == None:
    eps = mp.mpf(10)**(4 - mp.dps)
else:
    eps = mpmath.mpmathify(os.getenv('QUADRULE_EPS'))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

dim = int(sys.argv[1])
if dim < 1 or dim > 3:
    eprint('Invalid \'dim\' value!')
    exit(1)

# Read data from stdin
import subprocess
p = subprocess.Popen( \
        "sed -e 's/[[:space:]]\+/ /g' " + \
            "-e 's/^ *//g' -e 's/ *$//g' " + \
            "-e 's/ /,/g' " + \
            "-e 's/^/mpmath.mpmathify(\"/' -e 's/$/\")/' " + \
            "-e 's/,/\"),mpmath.mpmathify(\"/g' " + \
            "-e 's/\//\")\/mpmath.mpmathify(\"/g'",\
        shell=True, stdin=None, stdout=subprocess.PIPE)
qlist = []
n=0
for line in p.stdout.readlines():
    n += 1
    data = line.decode().replace('\n', '')
    if data == '':
        continue
    exec('point = [' + line.decode().replace('\n', '') + ']')
    if len(point) != dim + 1:
        eprint('Wrong number of columns at input data row', n)
        exit(1)
    qlist += [point]
if p.wait() != 0:
    eprint('Error processing input data, abort.')
    exit(1)

#-------- normalize the weights
total = mp.mpf(0)
for q in qlist:
    total += q[dim]
d = mp.mpf(1) - total
if d < -eps or d > eps:
    d = mp.mpf(1) / total
    for q in qlist:
        q[dim] *= d

#-------- optionally convert coordinates
if os.getenv('QUADRULE_VERTS') != None:
    s = os.getenv('QUADRULE_VERTS')
    fn = s.replace('file:', '')
    if fn != s:
        s = '`cat "'+fn+'"`'
    else:
        s = '"'+s+'"'
    p = subprocess.Popen('echo ' + s + ' | ' + \
            'sed -e "s/[][(){},;]/ /g" -e "s/  / /g"',
            shell=True, stdin=None, stdout=subprocess.PIPE)
    verts = []
    for line in p.stdout.readlines():
        for s in line.decode().replace('\n', '').split():
            verts += [mpmath.mpmathify(s)]
    if p.wait() != 0:
        eprint('Error occurred while processing QUADRULE_VERTS, abort.')
        exit(1)
    if len(verts) != dim * (dim + 1):
        eprint('Incorrect number of values (' + str(len(verts)) + \
                        ') defined by QUADRULE_VERTS, abort.')
        exit(1)
    # construct affine mapping: {v0,v_1,...v_d} -> {0,e_1,...,e_d}
    #       x -> A(x-v0)
    A = mp.inverse(mp.matrix([[verts[dim*(i+1)+j]-verts[j] \
                                for i in range(0,dim)] for j in range(0,dim)]))
    for q in qlist:
        x = A * mp.matrix([q[i] - verts[i] for i in range(0,dim)])
        for i in range(0,dim):
            q[i] = x[i,0]

#-------- output results, identifying and converting symmetry orbits

if os.getenv('QUADRULE_CSVOUT') != None and \
   os.getenv('QUADRULE_CSVOUT') != '0' and \
   os.getenv('QUADRULE_CSVOUT').upper() != 'N' and \
   os.getenv('QUADRULE_CSVOUT').upper() != 'NO':
    # output CSV format
    for q in qlist:
        tmp = str(q[0])
        for i in range(1, dim):
            tmp += ',' + str(q[i])
        print('Perm0(' + tmp + ')')
        print('Dup0(' + str(q[dim]) + ')')
    exit(0)

#--------- identify and convert symmetry orbits

if sys.version_info[0] >= 3:
    from functools import cmp_to_key    # for the 'cmp' keyword in sort()

def comp_float(a,b):
    c = a - b
    if c >= -eps and c <= eps:
        return 0
    elif c < 0:
        return -1
    else:
        return 1

def sorted_absc(p):
    q = p[0:dim]
    b = mp.mpf(1)
    for i in range(0, dim):
        b -= q[i]
    q += [b]
    if sys.version_info[0] >= 3:
        q.sort(key=cmp_to_key(comp_float))
    else:
        q.sort(cmp=comp_float)
    return q

def comp_list(n,m):
    # compare qlist[a] with qlist[b]
    p = sorted_absc(qlist[n])
    q = sorted_absc(qlist[m])
    for i in range(0, dim):
        k = comp_float(p[i], q[i])
        if k == 0:
            continue
        elif k < 0:
            return -1
        else:
            return 1
    return 0

if os.getenv('QUADRULE_P2O') != None and \
   os.getenv('QUADRULE_P2O') != '0' and \
   os.getenv('QUADRULE_P2O').upper() != 'N' and \
   os.getenv('QUADRULE_P2O').upper() != 'NO':
    # symmetry expansion: regard each point as a symmetry orbit
    if os.getenv('QUADRULE_P2O').upper() == "FIX_WEIGHTS":
        fix_weights = True
    else:
        fix_weights = False
    dups = []
    total = mp.mpf(0)   # sum of weights
    for q in qlist:
        p = sorted_absc(q)
        if dim == 1:
            if comp_float(p[0],mp.mpf(0.5)) == 0:
                print('Perm2(0.5)')
                dups += ['Dup2']
                total += q[dim]
            else:
                if fix_weights == True:
                    q[dim] /= mp.mpf(2)
                print('Perm11('+str(p[0])+')')
                dups += ['Dup11']
                total += mp.mpf(2) * q[dim]
        elif dim == 2:
            if comp_float(p[0],p[1]) == 0:
                if comp_float(p[0],mp.mpf(1)/mp.mpf(3)) == 0:
                    print('Perm3(1/3)')
                    dups += ['Dup3']
                    total += q[dim]
                else:
                    if fix_weights == True:
                        q[dim] /= mp.mpf(3)
                    print('Perm21('+str(p[0])+')')
                    dups += ['Dup21']
                    total += mp.mpf(3) * q[dim]
            elif comp_float(p[1],p[2]) == 0:
                if fix_weights == True:
                    q[dim] /= mp.mpf(3)
                print('Perm21('+str(p[1])+')')
                dups += ['Dup21']
                total += mp.mpf(3) * q[dim]
            else:
                if fix_weights == True:
                    q[dim] /= mp.mpf(6)
                print('Perm111('+str(p[0])+','+str(p[1])+')')
                dups += ['Dup111']
                total += mp.mpf(6) * q[dim]
        else:   # dim == 3
            if comp_float(p[0],p[1]) == 0 and comp_float(p[1],p[2]) == 0 and \
               comp_float(p[2],p[3]) == 0:
                print('Perm4(0.25)')
                dups += ['Dup4']
                total += q[dim]
            elif (comp_float(p[0],p[1]) == 0 and comp_float(p[1],p[2]) == 0) \
                 or \
                 (comp_float(p[1],p[2]) == 0 and comp_float(p[2],p[3]) == 0):
                if fix_weights == True:
                    q[dim] /= mp.mpf(4)
                print('Perm31('+str(p[1])+')')
                dups += ['Dup31']
                total += mp.mpf(4) * q[dim]
            elif comp_float(p[0],p[1]) == 0 and comp_float(p[2],p[3]) == 0:
                if fix_weights == True:
                    q[dim] /= mp.mpf(6)
                print('Perm22('+str(p[0])+')')
                dups += ['Dup22']
                total += mp.mpf(6) * q[dim]
            elif comp_float(p[0],p[1]) == 0:
                if fix_weights == True:
                    q[dim] /= mp.mpf(12)
                print('Perm211('+str(p[0])+','+str(p[2])+')')
                dups += ['Dup211']
                total += mp.mpf(12) * q[dim]
            elif comp_float(p[1],p[2]) == 0:
                if fix_weights == True:
                    q[dim] /= mp.mpf(12)
                print('Perm211('+str(p[1])+','+str(p[0])+')')
                dups += ['Dup211']
                total += mp.mpf(12) * q[dim]
            elif comp_float(p[2],p[3]) == 0:
                if fix_weights == True:
                    q[dim] /= mp.mpf(12)
                print('Perm211('+str(p[2])+','+str(p[0])+')')
                dups += ['Dup211']
                total += mp.mpf(12) * q[dim]
            else:
                if fix_weights == True:
                    q[dim] /= mp.mpf(24)
                print('Perm1111('+str(p[0])+','+str(p[1])+','+str(p[2])+')')
                dups += ['Dup1111']
                total += mp.mpf(24) * q[dim]
    # renormalize weights
    d = mp.mpf(1) / total
    for i in range(0, len(dups)):
        print(dups[i]+'('+str(qlist[i][dim]*d)+')')
    exit(0)

index = [*range(0,n)]
if sys.version_info[0] >= 3:
    index.sort(key=cmp_to_key(comp_list))
else:
    index.sort(cmp=comp_list)
i = 0
while i < n:
    j = i + 1
    while j < n and comp_list(index[i], index[j]) == 0 and \
        comp_float(qlist[index[i]][dim], qlist[index[j]][dim]) == 0:
        j += 1
    m = j - i
    #eprint("orbit length:", m)
    converted = False
    if dim == 1:
        if m == 1 and comp_float(qlist[index[i]][0], mp.mpf(0.5)) == 0:
            print('Perm2(0.5)')
            print('Dup2('+str(qlist[index[i]][dim])+')')
            converted = True
        elif m == 2 and comp_float(qlist[index[i]][0], mp.mpf(0.5)) != 0:
            print('Perm11('+str(qlist[index[i]][0])+')')
            print('Dup11('+str(qlist[index[i]][dim])+')')
            converted = True
    elif dim == 2:
        if m == 1 and comp_float(qlist[index[i]][0], mp.mpf(1)/mp.mpf(3)) == 0 \
                  and comp_float(qlist[index[i]][1], mp.mpf(1)/mp.mpf(3)) == 0:
            print('Perm3('+str(qlist[index[i]][0])+')')
            print('Dup3('+str(qlist[index[i]][dim])+')')
            converted = True
        elif m == 3:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[1]) == 0 or comp_float(p[1],p[2]) == 0:
                print('Perm21('+str(p[1])+')')
                print('Dup21('+str(qlist[index[i]][dim])+')')
                converted = True
        elif m == 6:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[1]) != 0 and comp_float(p[1],p[2]) != 0:
                print('Perm111('+str(p[0])+','+str(p[1])+')')
                print('Dup111('+str(qlist[index[i]][dim])+')')
                converted = True
    else:   # dim == 3
        if m == 1 and comp_float(qlist[index[i]][0], mp.mpf(0.25)) == 0 \
                  and comp_float(qlist[index[i]][1], mp.mpf(0.25)) == 0 \
                  and comp_float(qlist[index[i]][2], mp.mpf(0.25)) == 0:
            print('Perm4('+str(qlist[index[i]][0])+')')
            print('Dup4('+str(qlist[index[i]][dim])+')')
            converted = True
        elif m == 4:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[3]) != 0 and \
               ((comp_float(p[0],p[1]) == 0 and comp_float(p[1],p[2]) == 0) or \
                (comp_float(p[1],p[2]) == 0 and comp_float(p[2],p[3]) == 0)):
                print('Perm31('+str(p[1])+')')
                print('Dup31('+str(qlist[index[i]][dim])+')')
                converted = True
        elif m == 6:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[1]) == 0 and comp_float(p[2],p[3]) == 0 and \
                comp_float(p[0],p[2]) != 0:
                print('Perm22('+str(p[1])+')')
                print('Dup22('+str(qlist[index[i]][dim])+')')
                converted = True
        elif m == 12:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[1]) == 0 and \
               comp_float(p[1],p[2]) != 0 and comp_float(p[2],p[3]) != 0:
                print('Perm211('+str(p[1])+','+str(p[2])+')')
                print('Dup211('+str(qlist[index[i]][dim])+')')
                converted = True
            if comp_float(p[1],p[2]) == 0 and \
               comp_float(p[0],p[1]) != 0 and comp_float(p[2],p[3]) != 0:
                print('Perm211('+str(p[1])+','+str(p[0])+')')
                print('Dup211('+str(qlist[index[i]][dim])+')')
                converted = True
            if comp_float(p[2],p[3]) == 0 and \
               comp_float(p[0],p[1]) != 0 and comp_float(p[1],p[2]) != 0:
                print('Perm211('+str(p[2])+','+str(p[0])+')')
                print('Dup211('+str(qlist[index[i]][dim])+')')
                converted = True
        elif m == 24:
            p = sorted_absc(qlist[index[i]])
            if comp_float(p[0],p[1]) != 0 and comp_float(p[1],p[2]) != 0 and \
               comp_float(p[2],p[3]) != 0:
                print('Perm1111('+str(qlist[index[i]][0])+','\
                                 +str(qlist[index[i]][1])+','\
                                 +str(qlist[index[i]][2])+')')
                print('Dup1111('+str(qlist[index[i]][dim])+')')
                converted = True
    if not converted:
        for k in range(i,j):
            p = 'Perm'+str(dim+1)+'0('+str(qlist[index[k]][0])
            for l in range(1,dim):
                p += ','+str(qlist[index[k]][l])
            print(p+')')
            print('Dup0('+str(qlist[index[k]][dim])+')')
    i = j
