============
TIQC-SPICE
============

TIQC-SPICE/PyTIQC is a low-level physical simulation of Trapped Ion Quantum Computing (TIQC).

This project is inspired by SPICE, the famous simulator for integrated circuits developed at UC Berkeley in the 1970's. The goal is (eventually) to efficiently simulate the implementation of a quantum algorithm with a trapped ion quantum computer, and to be a useful tool for predicting system performance in the presence of noise and designing general quantum computers. 

TIQC-SPICE has its origins in Hartmut Haeffner's "quantumcomputer", a matlab suite for simulating TIQC experiments. Although the python code is re-written from scratch, the general code structure and many of the variable names remain.

Installation
------------

The following python packages are required. (Versions are the ones I've used for testing; higher versions will probably work, except that Parallel Python requires the server and client to have the same version). 

* Python 2.6.6
* Numpy 1.4.1
* Scipy 0.7.0
* Matplotlib 0.99.3
* Parallel Python 1.6.0
* python-cvxopt

After installing the above, add to path by (e.g.) adding this to your .bashrc::

  export PYTHONPATH="/home/youruserid/TIQC-SPICE/"

Files
-----

| PyTIQC/core         - main simulation engine
| PyTIQC/evaluation   - additional modules for analysis, usually in conjunction with data
| PyTIQC/experiments  - simulations of specific experiments / pulse sequences
| PyTIQC/testsuite    - unit test suite. also include simple examples of experiments
| PyTIQC/tools	    - additional tools, including external modules

How-To
------

The :doc:`tutorial` is a good place to start. Read this section for a step-by-step guide on writing and running a simulation.

The :doc:`libref` lists all the modules, classes, and functions defined in TIQC-SPICE. This portion of the documentation is automatically generated from docstrings.


Authors
-------

| Shannon X. Wang <sxwang@mit.edu>
| Thomas Monz <Thomas.Monz@uibk.ac.at>
| Richard Rines <rrines@mit.edu>

