#!/usr/bin/python
# File: sum_primes.py
# Author: VItalii Vanovschi
# Desc: This program demonstrates parallel computations with pp module
# It calculates the sum of prime numbers below a given integer in parallel
# Parallel Python Software: http://www.parallelpython.com

import math, sys, time
import pp

def isprime(n):
    """Returns True if n is prime and False otherwise"""
    if not isinstance(n, int):
        raise TypeError("argument passed to is_prime is not of 'int' type")
    if n < 2:
        return False
    if n == 2:
        return True
    max = int(math.ceil(math.sqrt(n)))
    i = 2
    while i <= max:
        if n % i == 0:
            return False
        i += 1
    return True

def sum_primes(n):
    """Calculates sum of all primes below given integer n"""
    return sum([x for x in xrange(2,n) if isprime(x)])

server_dict = {}
ppservers = []
for k in server_dict:
    ppservers.append(server_dict[k])
ppservers = tuple(ppservers)
if len(ppservers) > 0: ncpus = 0
else: ncpus = 'autodetect'

print ppservers

job_server = pp.Server(ncpus=ncpus, ppservers=ppservers, secret="tiqc_cluster1")

print "Starting pp with", job_server.get_ncpus(), "workers"

start_time = time.time()

# The following submits 8 jobs and then retrieves the results
inputs = [x for x in xrange(10000, 100000, 1000)]
jobs = [(input, job_server.submit(sum_primes,(input,), (isprime,), ("math",))) for input in inputs]
for input, job in jobs:
    print "Sum of primes below", input, "is", job()

print "Time elapsed: ", time.time() - start_time, "s"
job_server.print_stats()

print job_server.get_active_nodes()

# Parallel Python Software: http://www.parallelpython.com
