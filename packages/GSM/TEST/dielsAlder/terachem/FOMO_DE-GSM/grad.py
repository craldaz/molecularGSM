#!/global/software/Anaconda2/5.0.1/bin/python
import sys
import json
import numpy as np
import argparse
import os
#utility libraries
import pytc
import manage_xyz

# CHANGE THE CALC DETAILS HERE IF DOING SCF in method and grad_options
def dft(geom,atoms,fname,basis,isEnergy=False,isGradient=True):
    if (not (isinstance(basis, str))):
        raise TypeError(
            "basis must be a str")
    if (not (isinstance(isGradient, bool))):
        raise TypeError(
            "isGradient must be a bool")
    if (not (isinstance(isEnergy, bool))):
        raise TypeError(
            "isGradient must be a bool")
    #set the server up
    TC = pytc.Client(user='example', engine='Terachem', host='http://fire-05-31', verbose=False)
    TC.job_spec(atoms=atoms, charge=0, spinmult=1, closed_shell=True, restricted=True,  method='b3lyp', basis='6-31g')
    if isEnergy:
        results = TC.compute(name="E", job_type='energy', geoms=geom, **fomo_options)
        print(results.values()[0]['energy'][0])
    elif isGradient:
        grad_options=dict()
        results = TC.compute(name="Grad", job_type='gradient', geoms=geom, **grad_options)
        if isinstance(results.values()[0],dict):
            grad=results.values()[0]['gradient']
            energy=results.values()[0]['energy']
        else:
            count=0
            while not isinstance(results.values()[0],dict):
                results=TC.poll_for_results(results.keys()[0])
                count+=1
                if (count == 20):
                    print("JOB FAILED!!!")
                    return 0;
            grad=results.values()[0]['gradient']
            energy=results.values()[0]['energy']
    filename="scratch/GRAD"+fname+".1"
    with open(filename, 'w') as f:
        f.write(str(energy)+"\n")
        for N in range(len(atoms)):
            list=[grad[N*3+x] for x in range(3)]
            f.write(str(list)+"\n")

def fomo_cassci(geom,atoms,fname,closed,active,nstates,basis,wstate,readOrb,isEnergy=False,isGradient=False,isCoupling=False):
    # Check args are positive numbers and set their values
    if (not (isinstance(closed, int))) or (closed < 0):
        raise TypeError(
            "Must provide a positive int")
    if (not (isinstance(active, int))) or (active < 0):
        raise TypeError(
            "Must provide a positive int")
    if (not (isinstance(nstates, int))) or (nstates < 0):
        raise TypeError(
            "Must provide a positive int")
    if (not (isinstance(basis, str))):
        raise TypeError(
            "basis must be a str")
    if (not (isinstance(isGradient, bool))):
        raise TypeError(
            "isGradient must be a bool")
    if (not (isinstance(readOrb, bool))):
        raise TypeError(
            "readOrb must be a bool")
    if (not (isinstance(isCoupling, bool))):
        raise TypeError(
            "isGradient must be a bool")
    if (not (isinstance(isEnergy, bool))):
        raise TypeError(
            "isGradient must be a bool")
    #set the server up
    TC = pytc.Client(user='example', engine='Terachem', host='http://fire-05-31', verbose=False)
    TC.job_spec(atoms=atoms, charge=0, spinmult=1, closed_shell=True, restricted=True,  method='hf', basis='6-31g')

    # FOMO and CAS options
    fomo_options = {
        'casci':        'yes',
        'fon':          'yes',
        'closed':       closed,
        'active':       active,
        'cassinglets':  nstates,
        'precision': 'double',
        'threall': 1e-13,
        'convthre': 1e-7, 
        'purify': 'no',
        'maxit': 200,
    }
   
    orbitals="scratch/ORBFILE"+fname
    if (readOrb):
        try:
          with open(orbitals, "r") as orbfile:
              orbkey=orbfile.read()
        except IOError:
          print "Error: File does not appear to exist."
          return
        prev_calc=TC.poll_for_results(orbkey)
        prev_results=prev_calc.values()[0]
        orb_path = os.path.join(prev_results['job_scr_dir'], 'c0')

    if isEnergy:
        # Compute ground and excited state energies
        print fomo_options
        if readOrb:
            results = TC.compute(name="E", job_type='energy', geoms=geom, guess=orb_path,**fomo_options)
        else:
            results = TC.compute(name="E", job_type='energy', geoms=geom, **fomo_options)
        for i in range(nstates):
            print(results.values()[0]['energy'][i])

    elif isGradient and not isCoupling:
        grad_options = dict(fomo_options)
        results=[]
        grad_options['castargetmult'] = 1
        for n in wstate:
            grad_options['castarget']  = n-1
            if readOrb:
                results.append(TC.compute(name="Grad", job_type='gradient', geoms=geom,guess=orb_path, **grad_options))
            else:
                results.append(TC.compute(name="Grad", job_type='gradient', geoms=geom, **grad_options))
        #for n in range(len(wstate)):
        for n,w in enumerate(wstate):
            filename="scratch/GRAD"+fname+"."+str(w)
            if isinstance(results[n].values()[0],dict):
                grad=results[n].values()[0]['gradient']
                energy=results[n].values()[0]['energy']
            else:
                count=0
                while not isinstance(results[n].values()[0],dict):
                    results[n]=TC.poll_for_results(results[n].keys()[0])
                    count+=1
                    if (count == 10):
                        print("JOB FAILED!!!")
                        return 0;
                grad=results[n].values()[0]['gradient']
                energy=results[n].values()[0]['energy']
            with open(filename, 'w') as f:
                f.write(str(energy)+"\n")
                for N in range(len(atoms)):
                    list=[grad[N*3+x] for x in range(3)]
                    f.write(str(list)+"\n")
    elif isCoupling and not isGradient: 
        if (len(wstate) !=2):
            printf(" wstate must be 2")
            return 0
        nac_options = dict(fomo_options)
        nac_options['nacstate1'] = wstate[0]-1
        nac_options['nacstate2'] = wstate[1]-1
        nac_options['castargetmult'] = 1
        if readOrb:
            results = TC.compute(name="NAC", job_type='coupling', geoms=geom,guess=orb_path, **nac_options)
        else:
            results = TC.compute(name="NAC", job_type='coupling', geoms=geom, **nac_options)
        if isinstance(results.values()[0],dict):
            nac=results.values()[0]['nacme']
            energy=(results.values()[0]['energy'][wstate[0]-1] + results.values()[0]['energy'][wstate[1]-1])/2
        else:
            count=0
            while not isinstance(results.values()[0],dict):
                results=TC.poll_for_results(results.keys()[0])
                count+=1
                if (count == 10):
                    print("JOB FAILED!!!")
                    return 0;
            nac=results.values()[0]['nacme']
            energy=(results.values()[0]['energy'][wstate[0]-1] + results.values()[0]['energy'][wstate[1]-1])/2
        filename="scratch/COUP"+fname+"."+str(wstate[0]) + str(wstate[1])
        with open(filename, 'w') as f:
            f.write(str(energy)+"\n")
            for N in range(len(atoms)):
                list=[nac[N*3+x] for x in range(3)]
                f.write(str(list)+"\n")
    elif isGradient and isCoupling:
        grad_options = dict(fomo_options)
        results=[]
        #compute grad
        for n in wstate:
            grad_options['castarget'] = n-1
            grad_options['castargetmult'] = 1
            if readOrb:
                results.append(TC.compute(name="Grad", job_type='gradient', geoms=geom,guess=orb_path, **grad_options))
            else:
                results.append(TC.compute(name="Grad", job_type='gradient', geoms=geom, **grad_options))
        #compute nac
        nac_options = dict(fomo_options)
        nac_options['nacstate1'] = wstate[0]-1
        nac_options['nacstate2'] = wstate[1]-1
        nac_options['castargetmult'] = 1
        if readOrb:
            results.append(TC.compute(name="NAC", job_type='coupling', geoms=geom,guess=orb_path, **nac_options))
        else:
            results.append(TC.compute(name="NAC", job_type='coupling', geoms=geom, **nac_options))

        # get gradient
        for n in range(2):
            filename="scratch/GRAD"+fname+"."+str(n+1)
            if isinstance(results[n].values()[0],dict):
                grad=results[n].values()[0]['gradient']
                energy=results[n].values()[0]['energy']
            else:
                count=0
                while not isinstance(results[n].values()[0],dict):
                    results[n]=TC.poll_for_results(results[n].keys()[0])
                    count+=1
                    if (count == 20):
                        print("JOB FAILED!!!")
                        return 0;
                grad=results[n].values()[0]['gradient']
                energy=results[n].values()[0]['energy']
            with open(filename, 'w') as f:
                f.write(str(energy)+"\n")
                for N in range(len(atoms)):
                    list=[grad[N*3+x] for x in range(3)]
                    f.write(str(list)+"\n")
        # get coupling
        filename="scratch/COUP"+fname+"."+str(wstate[0]) + str(wstate[1])
        if isinstance(results[-1].values()[0],dict):
            nac=results[-1].values()[0]['nacme']
            energy=(results[-1].values()[0]['energy'][wstate[0]-1] + results[-1].values()[0]['energy'][wstate[1]-1])/2
        else:
            count=0
            while not isinstance(results[-1].values()[0],dict):
                results[-1]=TC.poll_for_results(results[-1].keys()[0])
                count+=1
                if (count == 10):
                    print("JOB FAILED!!!")
                    return 0;
            nac=results[-1].values()[0]['nacme']
            energy=(results[-1].values()[0]['energy'][wstate[0]-1] + results[-1].values()[0]['energy'][wstate[1]-1])/2
        with open(filename, 'w') as f:
            f.write(str(energy)+"\n")
            for N in range(len(atoms)):
                list=[nac[N*3+x] for x in range(3)]
                f.write(str(list)+"\n")
    try:
        orbkey = results[0].keys()[0].strip("[u'']")
        with open(orbitals, "w") as orbfile:
            orbfile.write(orbkey)
    except IOError:
        print 'error in writing orbs'
        return
    return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse Bool")   
    parser.add_argument('--fname', help='XYZ file', type=str, required=True)
    parser.add_argument('--method', help='FOMO or DFT', type=str, required=True)
    parser.add_argument('--active', help='number of active orbitals', type=int, required=False)
    parser.add_argument('--closed', help='closed orbital', type=int, required=False)
    parser.add_argument('--nstates', help='number of CI states', type=int, required=False)
    parser.add_argument('--wstate', type=int,nargs='+', help='<Required> Set flag', required=False)
    parser.add_argument('--basis', help='basis set', type=str, required=True)
    parser.add_argument("--readOrb", default=False, action="store_true" , help="Flag to read orb")
    parser.add_argument("--energy", default=False, action="store_true" , help="Flag to do something")
    parser.add_argument("--gradient", default=False, action="store_true" , help="Flag to do something")
    parser.add_argument("--coupling", default=False, action="store_true" , help="Flag to do something")
    args = parser.parse_args()
    
    #require active, closed, nstates else do not
    if args.method=="FOMO":
        if (args.active == None or args.closed == None or args.nstates==None):
            print( "ERROR!!")
            exit(1)
    geom1=manage_xyz.read_xyz("scratch/structure"+args.fname)
    atoms=manage_xyz.getAtoms(geom1)
    geom=manage_xyz.xyz_to_np(geom1)

    if args.method == "FOMO":
        fomo_cassci(geom,atoms,args.fname,args.closed,args.active,args.nstates,args.basis,args.wstate,args.readOrb,args.energy,args.gradient,args.coupling)
    if args.method == "DFT":
        dft(geom,atoms,args.fname,args.basis,isEnergy=False,isGradient=True)
