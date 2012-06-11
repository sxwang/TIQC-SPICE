=========
testsuite
=========

This directory contains a set of files for standard experiments which we often use in practice to characterize the system::

  testRabi.py
  testRamsey.py
  testSpectrum.py

The script ``testACStark.py`` runs checks the AC pulse the same way as for the Ramsey sequence. The script ``testACStark.py`` includes an optimization function to search for the value of the parameter ACcorr which gives the correct AC rabi time as observed in the experiment.

These scripts test some specific decoherence sources::

  testDecay.py
  testHeating.py

TODO: eventually the above should perhaps all be merged with:

testscripts
-----------

This is a set of unit tests using Python's unittest module. 

To use, run: testscripts/test_all.py.
