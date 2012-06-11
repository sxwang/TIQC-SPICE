==========
Parameters
==========

The class ``parameters`` contains variables that specifies the physical system and experimental parameters. All variables have a default value. These values are based on the Innsbruck Calcium-40 ion trap experiment. 

An instance of ``parameters`` is initialized by passing ``hspace`` as an argument::

    params = sim.parameters(hspace)

The standard ion trap parameters are defined in SI units as follows::

    params.omz = 2*pi*1.2  # secular frequency
    params.omc = 2*pi*0.023  # carrier freq in MHz
    params.omsb = 2*pi*0.023   # sideband freq in MHz
    params.pitime = pi / params.omc      # carrier pi-time
    params.omx = 2*pi*0.01  # corresponds to time for hiding/unhiding
    params.recoilangle_addr = 68  # addressed beam
    params.recoilangle_glob = 22 # global beam
    params.wavelength = 729e-9
    params.ionmass = 40 * 1.67e-27

Decoherence parameters are defined to the following standard values::

    # delay time between pulses due to beam switching/moving
    params.sameionswitchtime = 7
    params.diffionswitchtime = 46

    # spontaneous decay: lifetime of ca40 in us
    params.lifetime = 1168000

    # intensity fluctuation as a percentage
    params.intensityfluct = 0.02

    # coupling to spectator modes as a percentage
    params.specmodecoupling = 0.02

    # initialization error: probability of ending up in S'(m=+1/2)
    params.stateiniterr = 0.003

    # dephasing: times in us
    params.correlationTime = 333
    params.coherenceTime = 5000

    # heating: 1 quanta per "params.heatingrate" us
    params.heatingrate = 141000 

    # addressing error in percentage
    params.addressingerr = 0.035
    params.addressingerr_global = 0.02

    # hiding error 
    params.hidingerr = 0.99 # net population after each hide/unhide cycle
    params.hidingMeaserr = 0.005 # chance of ending up in D after each hide

Parallel python
~~~~~~~~~~~~~~~

``params`` expects a text file, ``server_dict`` in the core/ directory which contains a set of hostname-IP mappings of servers that can be used for parallelized cluster computing. The text file has the format::

  hostname: 192.168.5.1

lines beginning with ``#`` are interpreted as comments. If the file is not found, a warning is printed and the server list defaults to local only.

By default, only the local computer is used for parallelization and all of the cores are used. To select the actual servers to be used from all possible servers defined in ``server_dict``, use the function ``use_servers``::

    params.use_servers(['server0', 'server1']) #etc

If a server list is given, the local computer is not used. Alternatively, all available servers defined in ``server_dict`` as well as the local computer can be selected by::

    params.use_servers('all')

If using the local computer, the number of cores to be used can be set explicitly by::

    params.ppcpus = 'autodetect' # set to a number

By default, pp logs events to a file ``pp.log``. Logging can be disabled by setting::

    params.pplog = False

AC Stark pulse
~~~~~~~~~~~~~~

Here are the default values for the AC stark pulse. This pulse accomplishes a phase shift on a single ion by using both a detuning from the carrier transition and a residual stark shift due to the S-P transition. The detuning from the carrier is set to::

    params.ACdelta = -2*pi*20

The intensity of the laser is increased to compensate for the detuning::

    params.ACintenFac = 6.3

In the experiment, we measure the rabi frequency of the AC stark pulse to be::

    params.omac = 2*pi/80.*2

There is a correction factor that is set numerically to add the correct Stark shift factor to the Hamiltonian to account for the stark shift due to the S-P level, to obtain the correct rabi frequency::

    params.ACcorr = 0.0927

MS gate 
~~~~~~~~

All of the parameters related to the MS gate can be changed by the user. However, the method of setting up the parameters most closely matching what is done in the experiment is to select a value for the detuning. A positive value means detuned away from the carrier; a negative value means detuned toward the carrier. Unit is MHz::

    params.MSdelta = 2*pi*0.028

The parameter ``params.shortestMS`` refers to the denominator of the angle of the shortest MS gate in a pulse sequence. The default value is ``2``. For example, if a sequence contains the following pulses::

    sim.RMS(params, pi/2, 0)
    sim.RMS(params, pi/4, 0)
    sim.RMS(params, 3*pi/8, 0)

Then ``params.shortestMS`` should be set to 8. If the user does not explicitly set this value, it will be changed when a new MS pulse is instantiated.

The rest of the MS gate parameters are set to their optimal values by the function ``params.calcPulseParams``.

Use density matrix
~~~~~~~~~~~~~~~~~~

Since TIQC-SPICE is (so far) only meant to be used with states and not density matrices, an approximation need to be made in the case that we wish to simulate a pulse sequence with a density matrix as the initial state (for example, to simulate a partial sequence after some statistics have been obtained for a prior partial sequence). The density matrix input is decomposed into a list of eigenvectors. The eigenvectors with eigenvalues greater than 0.01 are collected into a dictionary and sequentially simulated. The results of these few simulations are summed with the eigenvalues as weights. Summing is done by the appropriate function in ``qctools``. In the parameter specification, a density matrix input can be specified by::

    params.use_rho0(rho0)

where ``rho0`` is a numpy array for the density matrix input.

Miscellaneous
~~~~~~~~~~~~~

A few other user-accessible parameters are listed here. 

An option is given to use the full Hamiltonian instead of the Hamiltonian with the Lamb-Dicke (LD) approximation. By default, LD approximation is turned on::

    params.LDapproximation = True

The default timestep for the unitary evolution is set to 1 us and the default timestep of the ODE solver is set to 0.01 us::

    params.stepsize = 1
    params.ODEtimestep = 0.01

During the evolution, the parameter ``printpulse`` control whether the pulse name, angle, phase, and ion are printed as each pulse is simulated; ``progbar`` controls whether the progressbar is displayed for each pulse. By default these are set to::

    params.printpulse = False
    params.progbar = True

The option to save data is set by::

    params.savedata = False

The default filename for the saved data is a timestamp followed by ``data-autosave.shlv``. Data is saved as python shelve objects.
