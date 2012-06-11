''' core of the simulation of a pulse sequence. '''

import numpy
import numpy as np
import scipy.linalg
import scipy.linalg as splg
import PyTIQC.tools.progressbar
import PyTIQC.tools.progressbar as progbar
import PyTIQC.tools.progressbar as progressbar

import simtools
import qmtools
import qctools
import pp
import sequel

import copy, logging, shelve

# these functions are only here for backward compatibility
def indexToState(ind, hspace):
    return qmtools.indexToState(ind,hspace)

def stateToIndex(statestr, hspace):
    return qmtools.stateToIndex(statestr, hspace)

# load and save entire runs
def loadRun(filename):
    """ load the data from a shelve object. Returns: pulseseq, params, dec, data """
    d = shelve.open(filename)
    data = d['data']
    params = d['params']
    pulseseq = d['pulseseq']
    dec = d['dec']

    return pulseseq, params, dec, data

def saveRun(pulseseq, params, dec, data, filename, clear=True):
    """ save all info from a simulation run into a shelve object. Remove the full time, state, and population arrays to save some space. """
    d = shelve.open(filename)
    if clear:
        if type(data) is list:
            for data_single in data:
                data_single.T = None
                data_single.Y = None
                data_single.YR = None
        else:
            data.T = None
            data.Y = None
            data.YR = None
    d['data'] = data
    d['params'] = params
    d['pulseseq'] = pulseseq
    d['dec'] = dec
    d.close()

#####################################
# core of the simulation
#####################################

def simulateevolution(pulseseq, params, dec):
    ''' a wrapper for simulateevolutionOnce, to run sim for varying starting states '''
    if dec.dict['all'] or dec.dict['switchtime']:
        pulseseq.addDelays(params) # add delays to account for switching

    if not params.y0_dict:
        data = simulateevolutionOnce(pulseseq, params, dec)
    else:
        data = None
        # if several starting states, iterate over them and sum at the end
        for [ev, y0] in params.y0_dict.items():
            params.y0 = y0
            data0 = simulateevolutionOnce(pulseseq, params, dec)

            if data == None:
                data = data0.times(ev)
            else:
                data0 = data0.times(ev)
                data = data.addstate(data0)

    return data


def simulateevolutionOnce(pulseseq, params, dec):
    """ a wrapper for simulationCore, to do randomizing internally
    by making changes to parameters and collecting the results """

    if dec.doSQL:
        sequel.insertSimToDB(pulseseq, params, dec)

    # the same time array as in simulationCore
    totaltime = pulseseq.totaltime
    T = np.append(np.arange(0, totaltime, params.stepsize), totaltime)

    if dec.doRandNtimes == 0:
        data = simulationCore(pulseseq, params, dec)
    else:
        k = 0

        ### addressing error
        if dec.dict['all'] or dec.dict['addressing']:
            params.set_addressing()
            for i in range(len(pulseseq)):
                pulseseq[i].targetion = params.addressing[pulseseq[i].ion]

        # save the variables we're going to change
        addressing_saved = np.copy(params.addressing)

        if dec.doPP:
            job_server = pp.Server( \
                ncpus = params.ppcpus, 
                ppservers = params.ppservers, 
                secret = params.ppsecret)

            if dec.params.pplog:
                ### job_server's logger
                job_server.logger.setLevel(logging.DEBUG)
                # create file handler which logs even debug messages
                fh = logging.FileHandler('pp.log')
                fh.setLevel(logging.DEBUG)
                # create formatter and add it to the handlers
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                fh.setFormatter(formatter)
                # add the handlers to logger
                job_server.logger.addHandler(fh)
                #######

            if dec.doPPprintstats:
                print "Starting PP with", job_server.get_ncpus(), \
                      "local workers and", params.ppservers

        if dec.progbar:
            widgets = [progbar.Percentage(), ' ', progbar.Bar(),' ', progbar.ETA() ]
            pbar = progbar.ProgressBar(widgets=widgets).start()

        # for now we store all param permutations, but be careful this could get very large!
        if dec.doPP:
            params_list = []
            pulseseq_list = []
            dec_list = []
            data = None
            pulseseq_orig = copy.deepcopy(pulseseq)

        while k < dec.doRandNtimes:

            #######################
            ### randomize the parameters

            ### initialization error
            if dec.dict['all'] or dec.dict['initerr']:
                rn = np.random.uniform()
                r_qubit = int(np.floor(np.random.uniform(0,params.hspace.nuions)))
                if rn < params.stateiniterr*params.hspace.nuions:
                    #print "   init error on ion ", r_qubit
                    params.addressing[:,r_qubit] = 0
                # propagate addressing matrix to the pulses
                for i in range(len(pulseseq)):
                    pulseseq[i].targetion = params.addressing[pulseseq[i].ion]

            ### spectator mode coupling, as initialized intensity shift
            if dec.dict['all'] or dec.dict['specmode']:
                rn = 1 + np.random.normal(scale=params.specmodecoupling)
                params.addressing = params.addressing*rn
                for i in range(len(pulseseq)):
                    pulseseq[i].targetion = params.addressing[pulseseq[i].ion]

            ### dephasing as offset in the pulse phase
            if dec.dict['phaseoffset']:
                phasenoise = np.random.normal(scale=params.phaseOffset, 
                                              size=len(pulseseq))
                for i in range(len(pulseseq)):
                    pulseseq[i].phase = pulseseq_orig[i].phase + phasenoise[i]

            ### dephasing
            if dec.dict['all'] or dec.dict['dephase']:
                dec.calcDephasing(T, params.stepsize)

            if len(pulseseq) != 0:
                ### spontaneous decay
                if dec.dict['all'] or dec.dict['spontdecay']:
                    stepsize = min(params.stepsize, pulseseq.totaltime/len(pulseseq))
                    dec.calcSpontaneousDecay(T, stepsize, dec.params.lifetime, params.hspace.nuions)

                ### heating
                if dec.dict['all'] or dec.dict['heating']:
                    stepsize = min(params.stepsize, pulseseq.totaltime/len(pulseseq))
                    dec.calcHeating(T, stepsize, dec.params.heatingrate)

                ### intensity fluct
                if dec.dict['all'] or dec.dict['intensity']:
                    stepsize = min(params.stepsize, pulseseq.totaltime/len(pulseseq))
                    dec.calcIntensFluct(T, stepsize, dec.params.intensityfluct)                

            #######################
            ### do it
            if not dec.doPP:
                data1 = simulationCore(pulseseq, params, dec)
            else:
                pulseseq_list.append(copy.deepcopy(pulseseq))
                params_list.append(copy.deepcopy(params))
                dec_list.append(copy.deepcopy(dec))

            #######################
            ### collect the results 

            # may have to adjust shape of result vector
            if not dec.doPP:
                if k == 0:
                    data = data1
                else:
                    try:
                        data += data1
                    except ValueError:
                        print "adding data failed, abandoning this instance"
                        continue

            ### initialization error
            # restore variables and increment counter
            if dec.dict['all'] or dec.dict['initerr'] or dec.dict['specmode']:
                params.addressing = np.copy(addressing_saved)
            k += 1

            ### update progressbar
            if dec.progbar and (not dec.doPP):
                pbar.update(int(1.*k*100/(dec.doRandNtimes)))

        if dec.doPP:   
            jobcounter = 0             
            runs = range(dec.doRandNtimes)
            for m in range( int(np.ceil( dec.doRandNtimes/float(dec.doRandBatch) )) ):
                if m < dec.doRandNtimes:
                    batch = runs[m*dec.doRandBatch : m*dec.doRandBatch+dec.doRandBatch]
                else:
                    batch = runs[m*doRandBatch:]

                jobs = [job_server.submit(simulationCore, \
                    args=(pulseseq_list[i], params_list[i], dec_list[i]), \
                    depfuncs=(), \
                    modules=('numpy','scipy', 'PyTIQC.core.simtools', \
                             'PyTIQC.core.qmtools', 'PyTIQC.core.sequel', \
                             'PyTIQC.tools.progressbar') ) \
                             for i in batch ]

                for job in jobs:
                    data1 = job()
                    if data1 == None:
                        print "simulationCore failed, continuing"
                        continue
                    if dec.progbar:
                        jobcounter+=1
                        pbar.update(int(1.*jobcounter*100/(dec.doRandNtimes)))
                    if not data:
                        data = data1
                    else:
                        try:
                            data += data1
                        except ValueError:
                            print "adding data failed, abandoning this instance"
                            continue
                
                if params.savedata:
                    saveRun(pulseseq, params, dec, data, params.savedataname, clear=False)

            # print pp stats
            if dec.doPPprintstats:
                print "PP server active: ", job_server.get_active_nodes()
                job_server.print_stats()

            #job_server.logger.removeHandler(fh)    
            if dec.params.pplog: logging.shutdown()
            job_server.destroy()
            
        # do averaging
        data.mean(k)

    return data

def simulationCore(pulseseq, params, dec):
    """ heart of the computation process. 
    input: pulse sequence, parameters (includes hilbert space def), and decoherence object. 
    output: a database object containing the time evolution of states. 
    """

    np = numpy
    splg = scipy.linalg
    pi = np.pi

    #qmtools = PyTIQC.core.qmtools
    #simtools = PyTIQC.core.simtools

    # make list of times and state vector
    totaltime = pulseseq.totaltime
    T = np.append(np.arange(0, totaltime, params.stepsize), totaltime)
    Y = np.zeros((len(T), len(params.y0)), np.complex128)
    Y[0,:] = params.y0

    # initialize indices and temp vars
    p0 = 0  # p0 is index to current pulse
    t0 = 0  # t0 is index to current time (in data list T)
            # tlen is amount of time to compute evolution
    ycur = np.mat(Y[0,:]).T  # convert to matrix to use *
    tcur = 0
    pcur = -1
    Ucur = 1

    # construct the time-dependent omrabi factors
    for pulse in pulseseq.seqqc:
        pulse.maketimedep(params.shape, params.doMSshapingCorr)

    # initialize the hamiltonian and noise objects
    ht = qmtools.Hamilton()
    ns = qmtools.Noise()
    

    # pre-fetch the noise dictionary
    if dec.doRandNtimes > 0:
        noise_dict = ns.Noise(params, dec)
        noise_total = [[noise_dict['none'][0]], [noise_dict['none'][1]], [noise_dict['none'][2]]]

        for key, [mult, add, uni] in noise_dict.iteritems():
            if (dec.dict['all'] or dec.dict[key]): # and not pulse.use_ideal:
                noise_total[0].append(mult)
                noise_total[1].append(add)
                noise_total[2].append(uni)

        noise_mult = ns.prodFunctions(noise_total[0])
        noise_add = ns.sumFunctions(noise_total[1])
        noise_uni = ns.sumFunctions(noise_total[2])

        projtimes = np.union1d(T[np.nonzero(dec.heatingV)], T[np.nonzero(dec.spontdecayV)])
    else:
        noise_mult = lambda t: 1
        noise_add = lambda t: params.hspace.operator_dict['zero']
        noise_uni = lambda t: params.hspace.operator_dict['zero']
        projtimes = np.array([])

    # storage of hidden ions
    hiddenions = np.zeros_like(params.addressing[-1])
    hiddenionsErr = np.ones_like(params.addressing[-1])
    hiddenionsCount = 0
    # a classical register to store intermediate measurements
    classical_reg = []

    if totaltime == 0:
        #set Y=y0 if no time for evolution
        Y[0,:] = params.y0
    else:
      ### time evolution starts here  
      while(tcur < totaltime):
        ### new pulse starts here  
        pulse = pulseseq.seqqc[p0]
        if pcur != p0:
            pcur = p0
            HTsaved = None
            if params.printpulse:
                print "Pulse ", p0, ":",  pulse
            if params.progbar and pulseseq.seqqc[p0].type != 'M':
                widgets = [progressbar.Percentage(), ' ', progressbar.Bar(),' ', progressbar.ETA()]
                pbar = progressbar.ProgressBar(widgets=widgets).start()
            #################
            # check for hiding. modify hiding matrix, then treat as delay
            nuions = params.hspace.nuions
            if pulse.type == "H":
                hiddenionsCount += 1
                if pulse.hide:
                    hiddenions[nuions-pulse.ion-1] = params.addressing[-1][pulse.ion]
                else:
                    hiddenions[nuions-pulse.ion-1] = 0
                    hiddenionsErr[nuions-pulse.ion-1] *= params.hidingerr
            else:
                pulse.targetion[np.nonzero(hiddenions)] = 0
                if dec.dict['all'] or dec.dict['hiding']:
                    # with some probability some ions are "lost" after unhiding
                    rn = np.random.uniform( size= \
                             len(np.nonzero(hiddenionsErr-1)[0]) )
                    for ith, ind in enumerate(np.nonzero(hiddenionsErr-1)[0]):
                        if rn[ith] > hiddenionsErr[np.nonzero(hiddenionsErr-1)[0][ith]]:
                            pulse.targetion[ind] = 0
                        else:
                            pulse.targetion = \
                                np.copy(params.addressing[pulse.ion])
                pulse.calculateIdealUnitary(params, np.nonzero(hiddenions[::-1])[0])
            # check for meas-init. measure and do projection, then treat as delay.
            if pulse.type == "I":
                if (dec.dict['all'] or dec.dict['hiding']) \
                        and pulse.incl_hidingerr:
                    reg = pulse.measure(np.asarray(ycur))
                    reg[1] += hiddenionsCount*params.hidingMeaserr
                    reg[0] -= hiddenionsCount*params.hidingMeaserr
                    classical_reg.append( reg )
                else:
                    classical_reg.append( pulse.measure(np.asarray(ycur)) )
                Uproj = pulse.Uinit()
                U = np.mat(Uproj)
                ynew = U * ycur
                ynew = ynew / np.sqrt(np.sum(np.power(abs(ynew),2)))
                ycur = ynew                

        #################
        # check for Uideal and skip time evolution if yes
        if pulse.use_ideal:
            ynew = np.dot(pulse.Uid, ycur)
            ycur = ynew
            tcur = pulse.endtime
            # save data and advance pulse
            t0 = np.searchsorted(T, tcur)
            if tcur == T[t0]:
                Y[t0,:] = np.asarray(ynew.T)
            else:
                T = np.insert(T, t0, tcur)
                [Yn1, Yn2] = np.array_split(Y, [t0])
                if len(Yn2) != 0:
                    Y = np.concatenate([Yn1, np.asarray(ynew.T), Yn2], axis=0)
                else:
                    Y = np.concatenate([Yn1, np.asarray(ynew.T)], axis=0)
            t0 = t0 + 1
            p0 = p0 + 1

        #################
        #### Unitary evolutions
        elif not pulse.dotimedepPulse and not params.dotimedep:

            assert tcur >= pulse.starttime \
                and tcur <= pulse.endtime \
                and abs(pulse.starttime + pulse.duration - pulse.endtime) < 0.001 \
                and pulse.duration >= 0, \
                "Pulse timing not consistent; missing copy.deepcopy?"

            tlen = min([pulseseq.seqqc[p0].endtime - tcur, T[t0+1]-tcur])

            [HT, lqc] = ht.Hamiltonian(pulse, params, LDApprox = params.LDapproximation)

            if not pulse.use_ideal:
                HT = noise_mult(tcur) * HT + noise_add(tcur)
                lqc = noise_mult(tcur) * lqc + np.diag(noise_add(tcur))


            Ugate = splg.expm2(-1j * tlen * HT)
            Ulqc1 = np.diag(np.exp(-1j * tcur * lqc))
            Ulqc2 = np.diag(np.exp(-1j * (-tlen-tcur) * lqc))

            U = np.mat(Ulqc2) * np.mat(Ugate) * np.mat(Ulqc1)
            ynew = U * ycur

            # normalize in case of jump operators for spontdecay and heating
            ynew = ynew / np.sqrt(np.sum(np.power(abs(ynew),2)))

            Ucur = U * Ucur

            # extra check for projective unitary (spontdecay, heating)
            if dec.doRandNtimes > 0 and np.sum(noise_uni(tcur)) != 0 \
                    and not pulse.use_ideal:
                U = np.mat(noise_uni(tcur))
                ynew = U * ycur
                ynew = ynew / np.sqrt(np.sum(np.power(abs(ynew),2)))
                ycur = ynew
                if abs(ynew[-1])**2 > 0.25:
                    print "Warning: further heating will exceed phonon space"

            ycur = ynew
            tcur = tcur+tlen

            # if reached data-point, store data and advance time
            datasaved = False
            if tcur == T[t0+1]:
                Y[t0+1,:] = np.asarray(ynew.T)
                t0 = t0 + 1
                datasaved = True

            if params.printpulse and params.progbar and pulseseq.seqqc[p0].type != 'M':
                pbar.update(int(1.*(tcur-pulseseq.seqqc[p0].starttime)*100 \
                                    /(pulseseq.seqqc[p0].duration)))

            # if pulse ended, advance to next pulse (if both then do both)
            if tcur == pulseseq.seqqc[p0].endtime:
                pulseseq.seqqc[p0].U = np.copy(np.asarray(Ucur))
                # save current data if not already saved
                if not datasaved:
                    T = np.insert(T, t0+1, tcur)
                    [Yn1, Yn2] = np.array_split(Y, [t0+1])
                    Y = np.concatenate([Yn1, np.asarray(ynew.T), Yn2], axis=0)
                    t0 = t0 + 1
                # advance to next pulse
                p0 = p0 + 1
                Ucur = 1

        #################
        #### Time-dependent HT: ODE solving
        else:
            # choose time step depending on detuning
            if pulseseq.seqqc[p0].detuning > 2*pi*2:
                stepduration = min(2/ (pulseseq.seqqc[p0].detuning/(2*pi)), 1)
            else:
                stepduration = params.ODEtimestep

            # for MS pulse, check and modify omrabi due to hiding
            if pulse.type == "M" and params.doMShidecorr:
                nuions = len(pulse.targetion)
                activeions = len(np.nonzero(pulse.targetion)[0])
                omc_fac = params.MShidecorr[activeions, nuions]
                if omc_fac == -1:
                    print "MS w/ hiding correction factor invalid, ignoring"
                else:
                    pulse.omrabi_b = params.omc_ms * omc_fac
                    pulse.omrabi_r = params.omc_ms * omc_fac

            # set up time-dep Hamiltonian
            if pulse.dobichro:
                HTblue = ht.Hamiltonian_timedep_complete(pulse.targetion, pulse.omrabi_bt, pulse.phase_light + pulse.phase_rb, pulse.detuning_b, params.omz, params.eta, params.hspace, LDApprox = params.LDapproximation)
                HTred =  ht.Hamiltonian_timedep_complete(pulse.targetion, pulse.omrabi_rt, pulse.phase_light - pulse.phase_rb, pulse.detuning_r, params.omz, params.eta, params.hspace, LDApprox = params.LDapproximation)
                HTorig = lambda t: HTblue(t) + HTred(t)
            else:
                HTorig = ht.Hamiltonian_timedep_complete(pulse.targetion, pulse.omrabi_t , pulse.phase, pulse.detuning, params.omz, params.eta, params.hspace, LDApprox = params.LDapproximation)
            HT = lambda t: noise_mult(t) * HTorig(t) + noise_add(t)

            psidot = lambda t,psi: -1j * np.dot(HT(t), psi)
            # ycur needs to be cast as a 1d array
            solver = qmtools.SEsolver(params.solver)

            Tloc = np.array([0.])
            Yloc = np.zeros([1,len(np.asarray(ycur))]) # this is to get the dimensions right, but make sure to remove the first row of 0's

            if params.progbar:
                widgets = [progressbar.Percentage(), ' ', progressbar.Bar(),' ', progressbar.ETA()]
                pbar = progressbar.ProgressBar(widgets=widgets).start()
            else:
                pbar = None

            # ODE solver
            # variable aliases for convenience
            pstarttime = pulseseq.seqqc[p0].starttime
            pendtime = pulseseq.seqqc[p0].endtime
            # first calculate the expected number of datapoints
            testtime = np.arange(tcur, pendtime, stepduration)
            testtime = np.append(np.delete(testtime, 0), pendtime)
            # extra check for projective unitary (spontdecay, heating)
            projtimes_cur = projtimes[np.intersect1d( \
                    np.nonzero(projtimes<pendtime)[0], \
                    np.nonzero(projtimes>pstarttime)[0]) ]
            for tproj in projtimes_cur:
                ms_ode = solver(HT, tcur, tproj, stepduration, np.ravel(ycur), \
                                pbar, pstarttime, pendtime)
                Tloc = np.append(Tloc, np.delete(ms_ode.time, 0))
                Yloc = np.append(Yloc, np.delete(ms_ode.y.transpose(), 0, axis=0), axis=0)
                tcur = tproj
                ycur = Yloc[-1,:].T
                if dec.doRandNtimes > 0 and np.sum(noise_uni(tcur)) != 0 \
                        and not pulse.use_ideal:
                    U = noise_uni(tcur)
                    ynew = np.dot(U, ycur)
                    ynew = ynew / np.sqrt(np.sum(np.power(abs(ynew),2)))
                    ycur = ynew
                    if abs(ynew[-1])**2 > 0.25:
                        print "Warning: further heating will exceed phonon space"
                    Yloc[-1,:] = np.asarray(ynew.T)
                
            ms_ode = solver(HT, tcur, pendtime, stepduration, np.ravel(ycur), \
                                pbar, pstarttime, pendtime)
            Tloc = np.append(Tloc, np.delete(ms_ode.time, 0))
            Yloc = np.append(Yloc, np.delete(ms_ode.y.transpose(), 0, axis=0), axis=0)

            Tloc = np.delete(Tloc, 0)
            Yloc = np.delete(Yloc, 0, axis=0)

            # check the length w.r.t. expected length of T vector and remove extras
            while len(Tloc) > len(testtime):
                for i in range(len(testtime)):
                    if testtime[i] != Tloc[i]:
                        Tloc = np.delete(Tloc, i)
                        Yloc = np.delete(Yloc, i, axis=0)
                        break

            if not params.saveallpoints:
                Tloc = [Tloc[-1]]
                Yloc = np.array([Yloc[-1]])

            # now we put the result into original result array
            # first, update the current time and state
            tcur = Tloc[-1]
            ycur = np.mat(Yloc[-1,:]).T
            # find the end of the pulse in the original T list
            t1 = np.nonzero(T >= tcur)[0][0]
            # only replace point if it's already been calculated
            if T[t1] == tcur:
                tend = t1+1
            else: # T[t1] > tcur
                tend = t1
            # replace the overlapping times in T with Tloc
            [Tnew1, Tnew2] = np.array_split( np.delete(T, range(t0+1, tend)) , [t0+1])
            Tnew = np.concatenate([Tnew1, Tloc, Tnew2])
            # mirror with Y and Yloc
            [Ynew1, Ynew2] = np.array_split( np.delete(Y, range(t0+1, tend), axis=0) , [t0+1])
            # seems that concatenate doesn't work on 2d arrays if they're empty
            if len(Ynew2) == 0:
                Ynew = np.concatenate([Ynew1, Yloc], axis=0)
            else:
                Ynew = np.concatenate([Ynew1, Yloc, Ynew2], axis=0)

            T = Tnew
            Y = Ynew

            p0 = p0 + 1
            t0 = np.nonzero(T >= tcur)[0][0]

    data = simtools.database(T,Y, params.hspace, pulseseq, register=classical_reg)

    data.creationtime = params.savedataname # timestamp the data

    if dec.doSQL:
        sequel.insertJobToDB(data)

    # get rid of lambda functions in order to send results back through pp
    for pulse in pulseseq.seqqc:
        pulse.omrabi_t = 0
        pulse.omrabi_bt = 0
        pulse.omrabi_rt = 0

    return data

# qctools.py ends here
