import re, os, sys, subprocess
import numpy
from pyscf import gto, scf
from pyscf.geomopt import geometric_solver, as_pyscf_method
from math import sqrt
import argparse

print('''
    Welcome to KirariNto. It is a program that tracks the electronic state during optimization, by comparing the NTO overlap between each 
    state with the target state of the previous geometry. It is based on a modified algorithm of Marco Campetella and Juan Sanz Garcia (10.1002/jcc.26162).
    Gaussian and ORCA are currently supported.
    A typical command:
    python3 kirariNto.py task.gjf
    Where task.gjf is an input file for TD-DFT single calculation of your initial geometry. nstates and nroot should be set explicitly.
    KirariNto invokes the program and Multiwfn, and the following environment variants are required:
    KIRARINTO_MULTIWFN
    KIRARINTO_GAUSSIAN
    KIRARINTO_FORMCHK
    KIRARINTO_ORCA
    KIRARINTO_ORCA2MKL
    Examples:
    export KIRARINTO_MULTIWFN=~/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn
    export KIRARINTO_GAUSSIAN=~/g16/g16
    export KIRARINTO_FORMCHK=~/g16/formchk
    export KIRARINTO_ORCA=/opt/orca5/orca
    export KIRARINTO_ORCA2MKL=/opt/orca5/orca_2mkl
    ''')

parser = argparse.ArgumentParser()
parser.add_argument('input', help='An input file for the TD-DFT calculation of your initial geometry.')
parser.add_argument('-c', '--constraints', help='(optional) An external geometry constraint file that follows the format requirement of PySCF.')
parser.add_argument('-p', '--program', help='The program you are using (Gaussian or ORCA; case insensitive).\
 Default: Gaussian if input file contains .gjf, and ORCA if input file contains .inp')
parser.add_argument('-d', '--delete', help='(optional) If the .chk or .gbw files to be deleted after each step. yes (default) or no.')
parser.add_argument('-g', '--grid', help='(optional) 1, 2 (default) or 3. The accuracy of integral grid for overlap evaluation. Denser for a large value.')
args = parser.parse_args()

INP = args.input
#PATH_MULTIWFN = '~/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn_noGUI'
PATH_MULTIWFN = '$KIRARINTO_MULTIWFN'
NUM_STEP = 1
GEOM = []
GEOM_OLD = []
LIST_ELEMENT = []
NROOT = 1
NSTATES = 3
HEADER = ''
TAIL = ''
NUM_ATOM = 0
PROG = ''
DELETE_CHK = True
GRID = 2

def popen(command):
    output, error = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).communicate()
    return output.split('\n')

def genNTO(fchkInput, logInput, indexOfState, fchkOutput):
    eigens = []
    f = popen(f'echo "18\n6\n{logInput}\n{indexOfState}\n2\n{fchkOutput}" | {PATH_MULTIWFN} {fchkInput}')
    isNTOEigens = False
    indexOfHOMO = 0
    for l in f:
        if re.match('\s*Orbitals from\s*\d+\s*to\s*\d+\s*are occupied',l):
            indexOfHOMO = int(l.split('to')[1].split('are')[0])
        if 'The highest 10 eigenvalues of NTO pairs' in l:
            isNTOEigens = True
        elif 'Sum of all eigenvalues' in l:
            isNTOEigens = False
        elif isNTOEigens:
            eigens.extend(l.split())
    eigens = [float(i) for i in eigens]
    return eigens, indexOfHOMO

def setOccForHole(fchkInput, moldenOutput, eigens, indexOfHOMO):
    currentOrbital = indexOfHOMO - 10
    commandStr = ''
    if currentOrbital > 0:
        commandStr = f'6\n26\n1-{currentOrbital}\n0\n'
    else:
        commandStr = f'6\n26\n'
    currentOrbital = indexOfHOMO
    for i in eigens:
        commandStr += f'{currentOrbital}\n{i ** 0.25 :.6f}\n'
        currentOrbital -= 1
    commandStr += f'00\n-1\n100\n2\n6\n{moldenOutput}\n'
    popen(f'echo "{commandStr}" | {PATH_MULTIWFN} {fchkInput}') # use read() to make the script wait for MULTIWFN

def setOccForElectron(fchkInput, moldenOutput, eigens, indexOfHOMO):
    commandStr = f'6\n26\n1-{indexOfHOMO}\n0\n'
    currentOrbital = indexOfHOMO + 1
    for i in eigens:
        commandStr += f'{currentOrbital}\n{i ** 0.25 :.6f}\n'
        currentOrbital += 1
    commandStr += f'00\n-1\n100\n2\n6\n{moldenOutput}\n'
    popen(f'echo "{commandStr}" | {PATH_MULTIWFN} {fchkInput}')

def CalcOverlap(molden1, molden2):
    f = popen(f'echo "5\n0\n1\n*,{molden2}\n1\n{GRID}" | {PATH_MULTIWFN} {molden1}')
    isIntegral = False
    for l in f:
        if 'Summing up all value and multiply differential element' in l:
            isIntegral = True
        elif isIntegral: 
            return float(l)

def compare(fchkA, fchkB, logA, logB, indeoxOfStateA, numberOfStates, prog = 'gaussian'):
    fileTail = '.fchk'
    if prog == 'orca':
        fileTail = '.molden'
    ntoFileA = fchkA.replace(fileTail, '_nto.fchk')
    holeFileA = fchkA.replace(fileTail, '_hole.molden')
    electronFileA = fchkA.replace(fileTail, '_electron.molden')
    eigensA, indexOfHOMO = genNTO(fchkA, logA, indeoxOfStateA, ntoFileA)
    setOccForHole(ntoFileA, holeFileA, eigensA, indexOfHOMO)
    setOccForElectron(ntoFileA, electronFileA, eigensA, indexOfHOMO)
    maxOverlap = 0
    maxOverlapIndex = 1
    for i in range(1, numberOfStates + 1):
        ntoFileB = fchkB.replace(fileTail, f'_nto_{i}.fchk')
        holeFileB = fchkB.replace(fileTail, f'_hole_{i}.molden')
        electronFileB = fchkB.replace(fileTail, f'_electron_{i}.molden')
        eigensB, indexOfHOMO = genNTO(fchkB, logB, i, ntoFileB)
        setOccForHole(ntoFileB, holeFileB, eigensB, indexOfHOMO)
        setOccForElectron(ntoFileB, electronFileB, eigensB, indexOfHOMO)
        overlapI = CalcOverlap(holeFileA, holeFileB) + CalcOverlap(electronFileA, electronFileB)
        print(f'NTO overlap between the state {i} and the previous TES is {overlapI:.5f}')
        if overlapI > maxOverlap:
            maxOverlap = overlapI
            maxOverlapIndex = i
    print(f'Now the state with largest overlap is determined to be {maxOverlapIndex}')
    for i in range(1, numberOfStates + 1):
        if i != maxOverlapIndex:
            popen(f"rm -rf {fchkB.replace('.fchk', f'_*_{i}.*')}")
    return maxOverlapIndex

def readForceAndGeomForGaussian(path):
    forceArr = []
    geomArr = []
    isNormalTermination = False
    with open(path) as f:
        isForce = False
        isGeom = False
        E = 0
        archivePart = ''
        isArchive = False
        for l in f.readlines():
            if 'Input orientation' in l:
                isGeom = True
                geomArr = []
            elif 'Distance matrix' in l or 'Rotational constants' in l:
                isGeom = False
            elif 'Forces (Hartrees/Bohr)' in l:
                isForce = True
            elif 'Cartesian Forces:  Max' in l:
                isForce = False
            elif 'SCF DONE' in l.upper():
                E = float(l.split('=')[1].upper().split('A.U.')[0])
            elif 'extrapolated energy' in l.lower():
                E = float(l.split('=')[1])
            elif 'E(TD-HF/TD-DFT)' in l.upper():
                E = float(l.split('=')[1])
            elif isArchive:
                # sometimes the MP2 energy may be separated by a \n. Two lines
                # have to be combined
                archivePart += l.upper().strip()
                isArchive = False
                E = float(archivePart.split('MP2=')[1].split('=')[0])
            elif isForce and (re.match("\\s*[0-9]+\\s+[0-9]+\\s*\\-*[0-9]+", l) is not None):
                forceArr.extend(l.split()[2:])
            elif isGeom and (re.match("\\s*[0-9]+\\s+[0-9]+\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                geomArr.extend(l.split()[3:])
            if 'Normal termination' in l:
                isNormalTermination = True
    if not isNormalTermination:
        raise Exception('Gaussian exited with error!')
    geomArr = [float(i) for i in geomArr]
    # forceArr = [float(i)/0.529 for i in forceArr]
    forceArr = [float(i) for i in forceArr]
    return [geomArr, forceArr, E]


def readForceAndGeomForORCA(path):
    forceArr = []
    geomArr = []
    splittedPath = path.split('.')
    # replace A.xxx into A.engrad
    path = ''.join(splittedPath[:-1]) + '.engrad'
    if os.path.exists(path): # for a vertical calculation, engrad file does not exist, and this will be jumped
        with open(path) as f:
            isForce = False
            isGeom = False
            isEnergy = False
            E = 0
            for l in f.readlines():
                if 'The atomic numbers and current coordinates in Bohr' in l:
                    isGeom = True
                    geomArr = []
                elif '#' in l and len(geomArr) > 3:
                    isGeom = False
                elif 'The current gradient' in l:
                    isForce = True
                elif '#' in l and len(forceArr) > 3:
                    isForce = False
                elif 'current total energy in Eh' in l:
                    isEnergy = True
                elif '#' in l and E != 0:
                    isEnergy = False
                elif isEnergy and (re.match("\\s*\\-*[0-9]+", l) is not None):
                    E = float(l.strip())
                elif isForce and (re.match("\\s*\\-*[0-9]+", l) is not None):
                    forceArr.append(l.strip())
                elif isGeom and (re.match("\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                    geomArr.extend(l.split()[1:])
    path = ''.join(splittedPath[:-1]) + '.out'
    with open(path) as f: 
    #the engrad file does not contain right TD-DFT energy for ORCA. Read it from .log
        for l in f.readlines():
            if 'E(tot)' in l:
                E = float(l.split('=')[1].strip().split()[0])
    geomArr = [float(i) * 0.52918 for i in geomArr]  # change Bohr to angstrom
    # ORCA outputs gradients. Here it is adapted to gaussian, which is the
    # force
    forceArr = [-float(i) for i in forceArr]
    return [geomArr, forceArr, E]

def readForceAndGeomForBAGEL(path, state):
    forceArr = []
    geomArr = []
    E = 0
    with open(path) as f:
        isForce = False
        isGeom = False
        isEnergy = False
        E = 0
        for l in f.readlines():
            if '*** Geometry ***' in l:
                isGeom = True
                geomArr = []
            elif 'Number of auxiliary basis functions' in l:
                isGeom = False
            elif 'Nuclear energy gradient' in l:
                isForce = True
            elif '* Gradient computed with' in l:
                isForce = False
            elif ' === FCI iteration ===' in l:
                isEnergy = True
            elif isEnergy and (re.match("\s*[0-9]+\s*[0-9]+\s*\*\s*\-*[0-9]+", l) is not None):
                if int(l.split()[1]) == int(state):
                    E = float(l.split()[-3])
            elif 'MS-CASPT2 energy : state' in l:
                if int(l.split('state')[1].split()[0]) == int(state):
                    E = float(l.split('state')[1].split()[1])
            elif isGeom and '{ "atom" :' in l:
                #{ "atom" : "C", "xyz" : [      7.990821,      1.210697,      3.574653 ] },
                thisAtom = l.split('[')[-1].split(']')[0].split(',')
                geomArr.extend(thisAtom)
            elif isForce and (re.match("\s*[x,y,z]\s*\-*[0-9]+", l) is not None):
                forceArr.append(l.split()[-1])
    # BAGEL outputs gradients. Here it is adapted to gaussian, which is the
    # force
    geomArr = [float(i) * 0.52918 for i in geomArr]  # change Bohr to angstrom
    forceArr = [-float(i) for i in forceArr]
    forceArr = addConstLag(geomArr, forceArr, CONSTRAINTS)
    if len(geomArr) == NUM_ATOM * 3:
        geomArr.extend(LAMBDAS)
    #print(f'Now reading energy for the state {state}, it is {E}')
    return [geomArr, forceArr, E]

def readForceAndGeom(path, state = 0):
    forceArr = []
    geomArr = []
    E = 0
    if PROG == 'gaussian':
        geomArr, forceArr, E = readForceAndGeomForGaussian(path)
    elif PROG == 'orca':
        geomArr, forceArr, E = readForceAndGeomForORCA(path)
    elif PROG == 'bagel':
        geomArr, forceArr, E= readForceAndGeomForBAGEL(path, state)
    else:
        raise Exception('Unsupported program!')
    return [geomArr, forceArr, E]

def parseGjfFile(path):
    geomArr = []
    elementArr = []
    header = ''
    tail = ''
    isGeom = False
    isHeader = True
    isTail = False
    with open(path) as f:
        for l in f.readlines():
            l = l.lower()
            if re.match("\s*\S+\s+\-*[0-9]*\.*[0-9]+\s+\-*[0-9]*\.*[0-9]+\s+\-*[0-9]*\.*[0-9]+\s*", l) is not None:
                if isHeader:
                    isHeader = False
                if (not isHeader) and (not isTail):
                    geomArr.extend(l.split()[1:])
                    elementArr.append(l.split()[0])
            elif isHeader:
                header += l
            elif not isHeader:
                isTail = True
                tail += l
    geomArr = [float(i) for i in geomArr]
    if (not 'nstates' in header) or (not 'root' in header):
        raise Exception('Please contain nstates and root options in your gjf file!')
    nstates = int(re.findall(r'nstates=(\d+)', header)[0])
    nroot = int(re.findall(r'root=(\d+)', header)[0])
    return [geomArr, elementArr, header, tail, nstates, nroot]

def parseORCAFile(path):
    geomArr = []
    elementArr = []
    header = ''
    tail = ''
    isGeom = False
    isHeader = True
    isTail = False
    with open(path) as f:
        for l in f.readlines():
            l = l.lower()
            if re.match("\s*\S+\s+\-*[0-9]*\.*[0-9]+\s+\-*[0-9]*\.*[0-9]+\s+\-*[0-9]*\.*[0-9]+\s*", l) is not None:
                if isHeader:
                    isHeader = False
                if (not isHeader) and (not isTail):
                    geomArr.extend(l.split()[1:])
                    elementArr.append(l.split()[0])
            elif isHeader:
                header += l
            elif not isHeader:
                isTail = True
                tail += l
    if (not 'nroots' in header) or (not 'iroot' in header):
        raise Exception('Please contain nroots and iroot options in the header part of your input file!')
    geomArr = [float(i) for i in geomArr]
    nstates = int(re.findall(r'nroots\s+(\d+)', header)[0])
    nroot = int(re.findall(r'iroot\s+(\d+)', header)[0])
    return [geomArr, elementArr, header, tail, nstates, nroot]

def writeGjf(Geom, Header, Tail, Name):
    f = open(Name, 'w+')
    f.write(Header)
    for i in range(NUM_ATOM):
        f.write('{ele}  {x:.8f}  {y:.8f}  {z:.8f}'.format(
            ele=LIST_ELEMENT[i], x=float(Geom[i * 3]), y=float(Geom[i * 3 + 1]), z=float(Geom[i * 3 + 2])))
        f.write('\n')
    f.write(Tail)
    f.write('\n')
    f.close()

def writeORCA(Geom, Header, Tail, Name):
    f = open(Name, 'w+')
    f.write(Header)
    f.write('\n')
    for i in range(NUM_ATOM):
        f.write('{ele}  {x}  {y}  {z}'.format(
            ele=LIST_ELEMENT[i], x=float(Geom[i * 3]), y=float(Geom[i * 3 + 1]), z=float(Geom[i * 3 + 2])))
        f.write('\n')
    f.write(Tail)
    f.write('\n')
    f.close()

def align(Geom1, Geom2): #transform Geom1 to align with Geom2
    xn = numpy.zeros([NUM_ATOM, 3])
    yn = numpy.zeros([NUM_ATOM, 3])
    geom1Aligned = Geom1[:]
    for i in range(NUM_ATOM):
        xn[i, 0:] =Geom1[i*3 : i*3+3]
        yn[i, 0:] =Geom2[i*3 : i*3+3]
    #Kabsch algorithm in 10.1107/S0567739476001873 for alignment
    R_matrix = yn.T @ xn
    RT_R = R_matrix.T @ R_matrix
    mius, A_matrix = numpy.linalg.eig(RT_R)
    mius = numpy.real(mius)
    mius = numpy.diag(mius)
    for i in range(3):
        mius[i, i] = 1/sqrt(mius[i, i])
    B_matrix = (mius @ (R_matrix @ A_matrix).T)
    U_matrix = B_matrix.T @ A_matrix
    xn = U_matrix @ xn.T
    for i in range(NUM_ATOM):
        geom1Aligned[i*3 : i*3+3] = xn.T[i, 0:]
    return geom1Aligned

def runEachStepForGaussian(mol):
    global NROOT
    global NUM_STEP
    global GEOM
    global GEOM_OLD

    GEOM_OLD = GEOM[:]
    GEOM = (mol.atom_coords() * 0.52917720859).flatten().tolist()
    geomAligned = []
    if NUM_STEP == 1:
        #For the first step, just add kwd "force" and run it
        thisHeader = HEADER.replace('#p', '#p force')
        geomAligned = GEOM[:]
    else: #otherwise, run a vertical calculation first
        thisHeader = HEADER.replace('#p', '#p guess=read')
        popen(f'cp {NUM_STEP - 1}.chk {NUM_STEP}.chk')
        geomAligned = align(GEOM, GEOM_OLD)

    thisHeader = thisHeader.replace('{chk}', f'{NUM_STEP}.chk')
    writeGjf(geomAligned, thisHeader, TAIL, f'JOBS/{NUM_STEP}.gjf')
    # at the vertical calculation step, first align the new molecule with the previous one
    popen(f'$KIRARINTO_GAUSSIAN JOBS/{NUM_STEP}.gjf; $KIRARINTO_FORMCHK {NUM_STEP}.chk')
    geomArr, forceArr, E = readForceAndGeom(f'JOBS/{NUM_STEP}.log')

    if NUM_STEP > 1:
        NROOT = compare(f'{NUM_STEP - 1}.fchk', f'{NUM_STEP}.fchk', f'JOBS/{NUM_STEP - 1}.log', f'JOBS/{NUM_STEP}.log', NROOT, NSTATES)
        thisHeader = HEADER.replace('#p', '#p force guess=read')
        thisHeader = re.sub(r'(root=)\d+', r'\g<1>' + str(NROOT), thisHeader)
        thisHeader = thisHeader.replace('{chk}', f'{NUM_STEP}.chk')
        # in the force calculation, use the unaligned geometry from PySCF
        writeGjf(GEOM, thisHeader, TAIL, f'JOBS/{NUM_STEP}.gjf')
        popen(f'$KIRARINTO_GAUSSIAN JOBS/{NUM_STEP}.gjf; $KIRARINTO_FORMCHK {NUM_STEP}.chk')
        geomArr, forceArr, E = readForceAndGeom(f'JOBS/{NUM_STEP}.log')
    popen(f'rm -rf {NUM_STEP - 1}_*.molden')
    popen(f'rm -rf {NUM_STEP - 1}_*.fchk')
    popen(f'rm -rf {NUM_STEP - 1}.fchk')
    if DELETE_CHK:
        popen(f'rm -rf {NUM_STEP - 1}.chk')
    NUM_STEP += 1
    return E, - numpy.array(forceArr).reshape(NUM_ATOM, 3)

def runEachStepForORCA(mol):
    global NROOT
    global NUM_STEP
    global GEOM
    global GEOM_OLD

    GEOM_OLD = GEOM[:]
    GEOM = (mol.atom_coords() * 0.52917720859).flatten().tolist()
    geomAligned = []
    if NUM_STEP == 1:
        #For the first step, just add kwd "engrad" and run it
        thisHeader = HEADER.replace('!', '! engrad')
        geomAligned = GEOM[:]
    else: #otherwise, run a vertical calculation first
        thisHeader = HEADER
        popen(f'cp {NUM_STEP - 1}.gbw {NUM_STEP}.gbw')
        geomAligned = align(GEOM, GEOM_OLD)

    writeORCA(geomAligned, thisHeader, TAIL, f'JOBS/{NUM_STEP}.inp')
    # at the vertical calculation step, first align the new molecule with the previous one
    popen(f'$KIRARINTO_ORCA JOBS/{NUM_STEP}.inp > JOBS/{NUM_STEP}.out; \
        cp JOBS/{NUM_STEP}.gbw .;\
        $KIRARINTO_ORCA2MKL {NUM_STEP} -molden;\
        mv {NUM_STEP}.molden.input {NUM_STEP}.molden')
    geomArr, forceArr, E = readForceAndGeomForORCA(f'JOBS/{NUM_STEP}.out')

    if NUM_STEP > 1:
        NROOT = compare(f'{NUM_STEP - 1}.molden', f'{NUM_STEP}.molden', f'JOBS/{NUM_STEP - 1}.out', f'JOBS/{NUM_STEP}.out', NROOT, NSTATES, prog='orca')
        thisHeader = HEADER.replace('!', '! engrad')
        thisHeader = re.sub(r'(nroots\s+)\d+', r'\g<1>' + str(NSTATES), thisHeader)
        thisHeader = re.sub(r'(iroot\s+)\d+', r'\g<1>' + str(NROOT), thisHeader)
        # in the force calculation, use the unaligned geometry from PySCF
        writeORCA(geomAligned, thisHeader, TAIL, f'JOBS/{NUM_STEP}.inp')
        popen(f'$KIRARINTO_ORCA JOBS/{NUM_STEP}.inp > JOBS/{NUM_STEP}.out; \
        cp JOBS/{NUM_STEP}.gbw .;\
        $KIRARINTO_ORCA2MKL {NUM_STEP} -molden;\
        mv {NUM_STEP}.molden.input {NUM_STEP}.molden')
        geomArr, forceArr, E = readForceAndGeomForORCA(f'JOBS/{NUM_STEP}.out')
    popen(f'rm -rf {NUM_STEP - 1}_*.molden')
    popen(f'rm -rf {NUM_STEP - 1}.molden')
    popen(f'rm -rf {NUM_STEP - 1}_*.fchk')
    if DELETE_CHK:
        popen(f'rm -rf {NUM_STEP - 1}.gbw')
    NUM_STEP += 1
    return E, - numpy.array(forceArr).reshape(NUM_ATOM, 3)

def runEachStep(mol):
    if PROG == 'gaussian':
        return runEachStepForGaussian(mol)
    elif PROG == 'orca':
        return runEachStepForORCA(mol)

def buildPySCFMolString():
    mol = ''
    for i in range(NUM_ATOM):
        mol += f'{LIST_ELEMENT[i]} {GEOM[i*3]} {GEOM[i*3+1]} {GEOM[i*3+2]};'
    return mol

def main():
    global GEOM
    global LIST_ELEMENT
    global NROOT
    global NSTATES
    global HEADER
    global TAIL
    global NUM_ATOM
    global PROG
    global DELETE_CHK
    global GRID

    if args.program != None:
        PROG = args.program.lower()
    else:
        if '.gjf' in INP.lower():
            PROG = 'gaussian'
        elif '.inp' in INP.lower():
            PROG = 'orca'
    if PROG == 'gaussian':
        GEOM, LIST_ELEMENT, HEADER, TAIL, NSTATES, NROOT = parseGjfFile(INP)
    elif PROG == 'orca':
        GEOM, LIST_ELEMENT, HEADER, TAIL, NSTATES, NROOT = parseORCAFile(INP)
    else:
        raise Exception('Unsupported program!')

    NUM_ATOM = int(len(GEOM) / 3)
    if '%chk' in HEADER:
        HEADER = re.sub(r'(%chk=).*?(\n)', r'\1{chk}\2', HEADER)
    elif PROG == 'gaussian':
        HEADER = '%chk={chk}\n' + HEADER
    if not '#p' in HEADER:
        HEADER = HEADER.replace('#', '#p') 
    popen('mkdir JOBS')
    
    mol = gto.M(buildPySCFMolString(), unit='Angstrom', verbose = 0)
    fake_method = as_pyscf_method(mol, runEachStep)

    params = {
    'convergence_energy': 1e-6,  # Eh
    'convergence_grms': 3e-4,    # Eh/Bohr
    'convergence_gmax': 4.5e-4,  # Eh/Bohr
    'convergence_drms': 1.2e-3,  # Angstrom
    'convergence_dmax': 1.8e-3  # Angstrom
    }
    if args.constraints != None:
        params['constraints'] = args.constraints
    if args.delete == 'no':
        DELETE_CHK = False
    if args.grid != None:
        GRID = int(args.grid)
    new_mol = fake_method.Gradients().optimizer(solver='geomeTRIC').kernel(params)


main()
#compare('u_td_1.fchk', 'u_td_2.fchk', 1,  3)
