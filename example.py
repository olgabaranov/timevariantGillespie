from dynamicGillespie import dynamicGillespie
import pickle
import networkx as nx
from numpy import logspace

R0=6
beta=1./10

gil = dynamicGillespie()  #initialise
gil.init_from_file('droplet.csv')   #read in links from file
gil.duration = 90   # set the duration of the simulation (in days)
gil.set_disease_params({'R0':R0,'beta':beta})
gil.simulate('original', int(gil.N*0.1), loop = 10) #run simulation 10 times (the loop argument)

path='results.pickle'
gil.save(path,mkdir = True) #save the results
