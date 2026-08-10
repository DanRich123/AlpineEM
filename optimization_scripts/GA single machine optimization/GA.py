from deap.tools.support import HallOfFame
import numpy as np
import random
from deap import creator,base, tools, algorithms
import json
import os
from matplotlib import pyplot as plt

#this is the def that runs your program/script of interest.
#It needs to take in values and output - return val,
def run_program(individual):
    #run fdtd
    params=[]
    for i in range(len(individual)):
        params.append(int(individual[i]))
    np.savetxt('params.txt', params)

    #run bias 1 and save data to a file
    os.system('python FDTD_full_integration_bias1.py')
    #evaluate the s parameters
    os.system('python Post_Processor_bias1.py')
    #load the s parameter data
    data=np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
    data=np.transpose(data)
    freq=data[0]
    s11_bias1_amp=data[1]
    s11_bias1_phase=data[2]

    #run bias 2 and save data to a file
    os.system('python FDTD_full_integration_bias2.py')
    #evaluate the s parameters
    os.system('python Post_Processor_bias2.py')
    #load the s parameter data
    data=np.loadtxt('S_parameters.csv', delimiter=',', skiprows=1)
    data=np.transpose(data)
    freq=data[0]
    s11_bias2_amp=data[1]
    s11_bias2_phase=data[2]

    #prep the goal we want
    #paremeters used in various ranges
    phase_resolution=15*np.pi/180
    s11_bias1_phase=np.unwrap(s11_bias1_phase)
    s11_bias2_phase=np.unwrap(s11_bias2_phase)
    phase_score=s11_bias1_phase-s11_bias2_phase

    #if statements are to avoid penalization for less/greater than desired goals
    #so here if you beat the 40 degree goal, it just returns a perfect score.
    #can add weights throughout to goals -if desired
    #least squares scoring
    count_11=0
    for i in range(len(freq)):
        if (freq[i]<=8 and freq[i]>=5):
            if phase_score[i]<phase_resolution:
                count_11+=(phase_score[i]-phase_resolution)**2/len(freq)
    #now score the s parameters based on whatever scoring criteria you want - i.e you can weight S11 or S21 totals if you want
    fitness_score=-1.0*(100*count_11)

    #plt.figure()
    #plt.plot(freq,phase_score)
    #plt.grid()
    #plt.savefig('current_phase_plot.png')

    #export the score
    return fitness_score,

#these are 2 custom def for me when I was using COMSOL with this script. These can be replaced with built in functions if desired.
#then I select two ranges of values that can be used.
def mu_uniform_seq(individual, sequence, indpb): 
    size=len(individual)
    for i, xseq in zip(range(size), sequence):
        if random.random() < indpb:
            individual[i] = random.choice(xseq)
    return individual,

def seq_func(start,stop,step):
    seq=[]
    counter=start-step
    while counter<(stop-step):
        counter=counter+step
        seq.append(round(counter,3))
    return seq

seq1=[4,5] #param of the material identifiers 1,2,3,4....

#Setup for GA
creator.create("FitnessMin", base.Fitness, weights=(1.0,)) #there's a multiweight option if desired.
creator.create("Individual", list, fitness=creator.FitnessMin)
toolbox=base.Toolbox()
toolbox.register("ep_1", random.choice, seq1)
toolbox.register("individual",tools.initRepeat,creator.Individual,(toolbox.ep_1),40) # tools.initCycle if >1 variable
toolbox.register("population", tools.initRepeat, list, toolbox.individual,20) #30 here is inital population size - how many individuals get run at once initially
toolbox.register("evaluate", run_program)
toolbox.register("mate", tools.cxUniform, indpb=0.15) #crossover probability for attributes within an individual
toolbox.register("mutate", mu_uniform_seq,sequence=[seq1], indpb=0.15) #mutate proability for attributes within an individual, use for num vars
toolbox.register("select", tools.selNSGA2)#tools.selTournament, tournsize=5) #selection method
#see https://deap.readthedocs.io/en/master/api/tools.html for operator options or see deap documation on github

#if desired, you can setup the initial population so it isn't generated randomly. Often referred to as seeding.
def initPopulation(pcls, ind_init, filename):
    contents = json.load(open(filename, "r"))
    return pcls(ind_init(c) for c in contents)
#toolbox.register("population_guess", initPopulation, list, creator.Individual, "my_guess.json")

# Initialize statistics object
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", np.mean, axis=0)
stats.register("std", np.std, axis=0)
stats.register("min", np.min, axis=0)
stats.register("max", np.max, axis=0)

#Initialize Algorithm variables
hof=tools.ParetoFront()

#use this is I want to use my seeds (initial guess)
#population=toolbox.population_guess()
#use this if no seeds (no initial guess), default option
population=toolbox.population()
#if you want to see what the initial population is before it runs
#print(population)

#this is the main loop, it will print as it goes
algorithms.eaMuPlusLambda(population, toolbox, mu=20,lambda_= 20, cxpb=0.7,mutpb=0.3,ngen=15,stats=stats,halloffame=hof,verbose=True)
#this prints the best parameters
#make a custom loop if we want it to print each gen while running.
print(hof[0])
np.save('best',hof)

#General Notes################################
#ngen give number of generation.
#population is the number of individuals per a generation.
#So total times your custom definition/program is run is (ngen)*population+initial population
#lambda is number of children to to produce at each geneartion - this is what is evaluated by the algorithm and produces a score
#mu is the numer of individuals to select for next generation population - pulls from both current population and current offspring in the case of eamupluslambda to produce the next round of lambda.
#this uses Varor but it's built in - it does mate or mutate or crossover each generation to produce lambda.
#cxpb is the probability that offspring is produced by crossover
#mutpb is the probability that offspring is produced by mutation
#These probabilities are different from 2 above in the setup section. These are whether or not the action takes places, and above is related to probablility that an attribute gets acted on within an individual.
#if we use MuCommaLambda we need to make Mu > Lambda, not sure how this works with NGSAII though if we don't use plus lambda.
