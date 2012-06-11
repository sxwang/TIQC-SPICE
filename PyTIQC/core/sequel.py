''' MySQL interface to jobs and simulation result records for TIQC-SPICE.
Needs to be implemented and incorporated into the main code. '''

import os, sys, time, string, re, glob, datetime
try:
    import MySQLdb, signal, datetime
    import PyTIQC.tools.tiqcSQL.models as sqlmod
except:
# if MySQL not installed then just ignore everything
    mysql=False
    pass

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'


def getsim(simid=None):
    # return a simulation set of parameters (name, param, pp)
    # given a simulation record ID (if None, then return first one)

    if not mysql:
        return {}

    dbx = MySQLdb.connect(user='tiqcspice',passwd='tiqc_cluster1',db='tiqcspice',host='marvin')
    db = dbx.cursor()
    
    sql = 'select * from Simulation'	
    db.execute(sql)
    for k in db.fetchall():
        (id, name, dt, param, pp, extra) = k
        if simid==None:
            break
        elif (id==simid):
            break
    db.close()
    dbx.close()
    print k
    return {'name':name, 'param':param, 'pp':pp, 'extra':extra}
    #records = sqlmod.Simulation.objects.all()
    #for rec in records:
    #    print rec


def insertSimToDB(pulseseq, params, dec):
    ''' create an entry for a Simulation '''

    if not mysql:
        return

    entry_ps = repr(pulseseq.seq)
    entry_params = MySQLdb.escape_string(repr(params.__dict__))
    entry_hspace = MySQLdb.escape_string(repr(params.hspace.__dict__))
    entry_dec = MySQLdb.escape_string(repr(dec.__dict__))

    dbx = MySQLdb.connect(user='tiqcspice',passwd='tiqc_cluster1',db='tiqcspice',host='marvin')
    db = dbx.cursor()
    	
    sql = "insert into Simulation (name, pulseseq, params, hspace, decoherence) values ('%s', '%s','%s','%s','%s')" \
        % (dec.doSQLname, entry_ps, entry_params, entry_hspace, entry_dec)
    try:
        db.execute(sql)
    except Exception, e:
        print "ERROR in sql insertSimToDB:", e

    db.close()
    dbx.close()


def insertJobToDB(data):
    ''' unpack the database object and insert into DB '''
    if not mysql:
        return

    entry_T = repr(data.TP)
    entry_Y = repr(data.YP)
    entry_hspace = MySQLdb.escape_string(repr(data.hspace.__dict__))

    dbx = MySQLdb.connect(user='tiqcspice',passwd='tiqc_cluster1',db='tiqcspice',host='marvin')
    db = dbx.cursor()
    	
    sql = "insert into Job (time, state, hspace) values ('%s','%s','%s')" \
        % (entry_T, entry_Y, entry_hspace)
    try:
        db.execute(sql)
    except Exception, e:
        print "ERROR in sql insertJobToDB:", e

    db.close()
    dbx.close()    
    

if __name__=="__main__":
    print "getting a simulation specification"
    getsim()

