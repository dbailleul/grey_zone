############################################################################
####Modèle 3 populations et Bottleneck, test
############################################################################

import simuPOP as sim

from simuPOP.sampling import drawRandomSample
from simuPOP.demography import *

import math
import numpy as np
import os
import time

global NOM
NOM = time.strftime('%Y-%m-%d-%Hh %Mmin',time.localtime())
print NOM

global CHEMIN
CHEMIN = os.getcwd()
print CHEMIN

def reccord (chaine,nfic):
    cheminfst = os.path.join(CHEMIN, (NOM + nfic + '.txt'))
    resultsfst=open(cheminfst,'a')
    resultsfst.write(chaine)
    resultsfst.close
    return cheminfst

def calcFst(pop):
    sortie = ''
    sim.stat(pop, structure=range(10), vars=['F_st'])
    Fstpop = pop.dvars().F_st
    for a in range(100):
        sample = drawRandomSample(pop, sizes=[50]*2)
        sim.stat(sample, structure=range(10), vars=['F_st'])
        Fstsample = sample.dvars().F_st
        sample.addInfoFields('order')
        order = list(range(100))
        fstsim = ''
        for rep in range(1000):
            merged=sample
            merged.mergeSubPops()
            np.random.shuffle(order)
            merged.setIndInfo(order, field = 'order')
            merged.sortIndividuals('order')
            merged.splitSubPop(0, [50]*2)
            sim.stat(merged, structure=range(10), vars=['F_st'])
            fstsim += '%s\t' % merged.dvars().F_st
        sortie += '%3d\t%.6f\t%3d\t%.6f\t%s\n' % (pop.dvars().gen, Fstpop, a, Fstsample, fstsim)
    reccord (sortie, "dataout")
    return True

model = EventBasedModel(
...     N0=([1000000]*3),
...     events=[
...         ResizeEvent(at=500, sizes=1000),
...     ]
... )

pop = sim.Population(size=model.init_size, loci=[10]*10,
...     infoFields='migrate_to')

pop.evolve(
initOps=[
sim.InitSex(),
sim.InitGenotype(freq=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1])
],
preOps=[
	sim.Migrator(rate=[
		[0,0.000005,0.000005],
		[0.000005,0,0.000005],
		[0.000005,0.000005,0],
	],mode=sim.BY_PROPORTION, end=499),
	sim.Migrator(rate=[
		[0,0.005,0.005],
		[0.005,0,0.005],
		[0.005,0.005,0],
	],mode=sim.BY_PROPORTION, begin=500)
],
matingScheme=sim.RandomMating(subPopSize=model,ops=sim.Recombinator(rates=0.01)),
postOps=[
sim.PyOperator(func=calcFst, step=100)
],
gen = 1000
)

for idx, name in enumerate(pop.subPopNames()):
...     print('%s (%d): %.4f' % (name, pop.subPopSize(name),
...         pop.dvars(idx).alleleFreq[0][0]))
