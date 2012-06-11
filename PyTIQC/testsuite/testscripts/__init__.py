
#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-30 23:03:15 tomi"

#  file       __init__.py
#  author     Thomas Monz 2011

"""A testing framework for the TIQC-SPICE/PyTIQC package
"""
import unittest


#------------------------------------------------------------------------------
class Test_BaseFunctions(unittest.TestCase):
    pass
#------------------------------------------------------------------------------
# Collect all test suites for running
all_suites = unittest.TestSuite((
  unittest.makeSuite(Test_BaseFunctions)
  ))

plot_suites = unittest.TestSuite((
  unittest.makeSuite(Test_BaseFunctions)
  ))

# Run all sub-test modules in this package by importing them
import test_carrier
test_carrier.TestUserFunction.previous_tests = all_suites.countTestCases()
all_suites.addTest(test_carrier.all_suites)
plot_suites.addTest(test_carrier.plot_suites)

import test_evaluations
test_evaluations.TestUserFunction.previous_tests = all_suites.countTestCases()
all_suites.addTest(test_evaluations.all_suites)
plot_suites.addTest(test_evaluations.plot_suites)

import test_decoherence
test_decoherence.TestUserFunction.previous_tests = all_suites.countTestCases()
all_suites.addTest(test_decoherence.all_suites)
plot_suites.addTest(test_decoherence.plot_suites)

def run_all():
    unittest.TextTestRunner(verbosity=2).run(all_suites)

def run_plot():
    test_carrier.TestUserFunction.figNum = 0
    test_evaluations.TestUserFunction.figNum = 0
    test_decoherence.TestUserFunction.figNum = 0
    unittest.TextTestRunner(verbosity=2).run(plot_suites)

def debug():
    all_suites.debug()


# __init__.py ends here
