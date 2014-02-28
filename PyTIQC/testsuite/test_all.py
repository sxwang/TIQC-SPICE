#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-27 14:44:54 c704252"

#  file       test_all.py
#  author     Thomas Monz 2011

import testscripts
import sys

if __name__ == "__main__":
#    testscripts.run_all()
    try:
        if sys.argv[1] == "plot":
            testscripts.run_plot()
    except IndexError:
        testscripts.run_all()


# test_all.py ends here
