class dynamicGillespie():

    def __init__(self):

        self.linkdict = None
        self.N = 0
        self.G = None
        self.users = None
        self.user_idx = None
        self.stepsize = None
        self.duration_data = None
        self.duration = None
        self.state = None
        self.neighbors = None
        self.beta = None
        self.alpha = None
        self.R0 = None
        self.s_recovery = {}
        self.s_infection = {}
        self.s_lambda = {}
        self.s_linkweight = None
        self.results = []

    def init_from_file(self, path_to_file, daynight = False):
        """takes a csv file as input; format [user1 user2 timestamp duration], will not work if the duration is not equal for all links"""
        links = self.read_file(path_to_file)
        self.init_from_list(links, daynight = daynight)

    def init_from_list(self, linklist, daynight = False):
        """takes a list of links as input; format [user1 user2 timestamp duration], will not work if the duration is not equal for all links"""
        from data_processing import list_to_graph
        if daynight: from data_processing import list_to_daynight

        if len(set([x[3] for x in linklist]))>1:
            print "incompatible list format; time slices must have the same duration throughout the list"
            return

        self.users = list(set([y for x in linklist for y in [x[0],x[1]]]))
        self.user_idx = [idx for idx, x in enumerate(self.users)]
        self.N = float(len(self.users))
        self.stepsize = float(linklist[0][3]) / (60 * 60 * 24.)
        self.duration_data = float(max([x[2] for x in linklist]) + self.stepsize - min([x[2] for x in linklist])) / (60 * 60 * 24.)
        self.duration = 90

        self.linklist_to_dict(linklist)
        self.G = list_to_graph(linklist, self.duration_data, self.stepsize, nodes_as_int = False)
        if daynight:
            from datetime import datetime as dt
            Glist = list_to_daynight(linklist,[8,20],nodes_as_int = False, exportIDdict=False)
            self.G_day = Glist[0]
            self.G_night = Glist[1]
            self.starting_epsilon = 1 - dt.fromtimestamp(float(sorted(linklist,key=lambda x: x[2])[0][2])).hour/24.

        self.neighbors = {x : [] for x in self.users}
        self.state = {x : 0 for x in self.users}
        self.s_recovery = {x : 0 for x in self.users}
        self.s_infection = {x : 0 for x in self.users}

    def read_file(self, path_to_file):
        """ invoked by init_from_file"""
        #read file
        with open(path_to_file, 'r') as f:
            link_list = [line.rstrip('\r\n').split(' ') for line in f]
        # remove selflinks
        links=[[x[0],x[1],int(x[2]),int(x[3])] for x in link_list[1:] if not x[0]==x[1]]

        return links

    def linklist_to_dict(self, linklist):
        keys = set([x[2]/ (60 * 60 * 24.) for x in linklist])
        tmin = min(keys)
        self.linkdict = {x-tmin:[] for x in keys}
        for edge in linklist:
            self.linkdict[edge[2] / (60 * 60 * 24.) - tmin].append((edge[0],edge[1]))


    def calculate_rho(self,timestep = None):
        """aka overall mean link weight"""
        if timestep:
            links = self.linkdict[timestep]
            rho = (2 * len(links)) / (self.N * (self.N - 1))
        else:
            rho = self.calc_meandegree() / float(self.N-1)
        return rho

    def calc_alpha(self,R0,beta):
        return R0*beta/self.calc_meandegree()

    def calc_beta(self,R0,alpha):
        return alpha*self.calc_meandegree()/R0

    def calc_R0(self,alpha,beta):
        return (alpha/beta)*self.calc_meandegree()

    def calc_meandegree(self):
        """mean degree of the time averaged network"""
        import networkx as nx
        return sum([x[2]['weight'] for x in nx.DiGraph(self.G).edges(data=True)])/float(len(self.G.nodes()))

    def set_disease_params(self, param_dict):
        """functions which makes sure that R0, alpha and beta are consistent"""
        keyset = set(param_dict.keys())
        if len(keyset) < 2:
            print "provide at least two of the three [alpha, beta, R0]"
        elif set(["R0", "alpha", "beta"]).issubset(keyset):
            assert param_dict["R0"]==self.calc_R0(param_dict["alpha"],param_dict["beta"])
        elif set(["R0", "beta"]).issubset(keyset):
            param_dict["alpha"] = self.calc_alpha(param_dict["R0"],param_dict["beta"])
        elif set(["R0", "alpha"]).issubset(keyset):
            param_dict["beta"] = self.calc_beta(param_dict["R0"],param_dict["alpha"])
        elif set(["alpha", "beta"]).issubset(keyset):
            param_dict["R0"] = self.calc_R0(param_dict["alpha"],param_dict["beta"])

        self.R0 = param_dict["R0"]
        self.alpha = param_dict["alpha"]
        self.beta = param_dict["beta"]

    def reset_system(self):
        """reset all the internal stuff back to initial state"""
        self.state = {user : 0 for user in self.users}
        self.neighbors = None
        self.s_recovery = {user : 0 for user in self.users}
        self.s_infection = {user : 0 for user in self.users}
        self.s_lambda = {user : 0 for user in self.users}


    def infect_random_users(self,victims):
        """victims can be either an int, then the infected are picked at random, or a list with node ids"""
        if type(victims) == float: victims = int(victims)
        if type(victims) == int:
            from random import sample
            infected = sample(self.state.keys(),victims)
            for user in infected:
                self.state[user] = 1
        elif type(victims) == list:
            for user in victims:
                self.state[user] = 1

    def update_neighbors(self, timestep, how):
        """invoked each timestep to update the dictionary of neighbors; how is the simulation mode"""
        if how == 'meanfield':
            self.meanfield_weights()
        elif how == 'original':
            self.original_weights(timestep)
        elif how == 'meanfield_modulated':
            self.s_linkweight = self.calculate_rho(timestep)
            self.meanfield_weights()
        elif how == 'the_other':
            self.the_other_weights(timestep)
        elif how == 'averaged':
            self.averaged_weights()
        elif how == 'daynight':
            self.daynight_weights(timestep)

    def meanfield_weights(self):
        """connects every two nodes by a link with time-link-averaged value for the network"""
        self.neighbors = {user : {} for user in self.users}
        for user1 in self.users:
            for user2 in self.users:
                if not user1 == user2:
                    self.neighbors[user1][user2] = self.s_linkweight
                    self.neighbors[user2][user1] = self.s_linkweight

    def averaged_weights(self):
        """takes the link weights from the time averaged network"""
        self.neighbors = {user : {} for user in self.users}
        for edge in self.G.edges(data=True):
            self.neighbors[edge[1]][edge[0]] = edge[2]['weight']

    def daynight_weights(self,time):
        """the weights are pulled from either day (time averaged 8 a.m to 8 p.m.) or night (time averaged 8 p.m. to 8 a.m.) activity networks"""
        if time%1 >= 8./24 and time%1 < 20./24:
            G = self.G_day
        else:
            G = self.G_night
        self.neighbors = {user : {} for user in self.users}
        for edge in G.edges(data=True):
            self.neighbors[edge[1]][edge[0]] = edge[2]['weight']

    def original_weights(self, timestep):
        """takes links according to the link list, weight is always 1"""
        links = self.linkdict[timestep]
        self.neighbors = {user : {} for user in self.users}
        for edge in links:
            self.neighbors[edge[0]][edge[1]] = self.s_linkweight
            self.neighbors[edge[1]][edge[0]] = self.s_linkweight

    def the_other_weights(self, timestep):
        """ mean field time-link-averaged links modulated by the current activity of the network (the purple approximation) """
        self.neighbors = {user : {} for user in self.users}
        rho_over_rho = self.calculate_rho(timestep) / self.rho
        for edge in self.G.edges(data=True):
            self.neighbors[edge[1]][edge[0]] = edge[2]["weight"] * rho_over_rho

    def calculate_recovery_process(self, user):

        #assigns the recovery rate according to the state which is in the dict ATM!
        if self.state[user]:
            self.s_recovery[user] = self.beta
        else:
            self.s_recovery[user] = 0

    def calculate_infection_process(self):
        for user in self.users:
            if self.state[user]:
                self.s_infection[user] = 0
            else:
                self.s_infection[user] = sum([weight * self.state[nbr] for nbr,weight in self.neighbors[user].items()]) * self.alpha

    def update_lambda(self):
        self.s_lambda = {user : self.s_recovery[user] + self.s_infection[user] for user in self.users}

    def update_state(self, node, time):
        """ changes the state of the network """
        if self.state[node]:
            self.state[node] = 0
        else:
            self.state[node] = 1

    def update_results(self, node, ISlinks, time):
        """ appends the current prevalence and time in the results vector"""
        #this function have to be invoked AFTER! the state of the node was updated
        self.results[-1]['timeline'].append(time)

        if self.state[node]: self.results[-1]['prevalence'].append(self.results[-1]['prevalence'][-1] + 1)
        else: self.results[-1]['prevalence'].append(self.results[-1]['prevalence'][-1] - 1)
        self.results[-1]['ISlinks'].append(ISlinks)

    def check_before_simulation(self):
        """ atm redundand """
        return all([self.users,
                    self.linkdict,
                    self.stepsize,
                    self.duration,
                    # type(self.beta) == int or type(self.beta)== float,
                    # type(self.alpha) == int or type(self.alpha)== float,
                    self.state,
                    self.s_recovery,
                    self.s_infection])

    def countISLinks(self):
        #calculates the number of IS  links in the current state
        infs = [key for key, value in self.state.items() if value ==1]
        ISlinks = 0
        #import pdb; pdb.set_trace()
        ISlinks += sum([weight for inf in infs for nbr,weight in self.neighbors[inf].items() if self.state[inf]+self.state[nbr] == 1 ])
        return ISlinks

    def simulate(self, mode, infection_seed, loop = 1):
        """ the main simulation function; infection seed can be eiter an int (then the seeds are chosen randomly) or a list with node IDs; modes can be chosen from: original, meanfield, meanfield_modulated, averaged, daynight, the_other """
        from scipy.stats import rv_discrete
        import networkx as nx
        import numpy as np

        if not self.check_before_simulation():
            print "check if the initialisation was complete"
            return

        #from numpy.random import exponential

        if mode == 'meanfield':
            self.s_linkweight = self.calculate_rho()
        elif mode == 'original':
            self.s_linkweight = 1 # check if this is really the case!!!!
        elif mode == 'the_other':
            self.rho = self.calculate_rho()
        if mode == 'daynight':
            self.regular_stepsize = self.stepsize
            self.stepsize = 0.5


        counter = 0
        while counter < loop:

            #draw initial tau
            tau = np.random.exponential(1)
            epsilon = 1
            cur_time = 0

            self.reset_system()
            self.infect_random_users(infection_seed)
            #nessesary here to count the IS links at the beginning
            if mode == 'daynight':
                self.update_neighbors(cur_time, mode)
            else:
                timestep = sorted([(x,abs(cur_time % self.duration_data - x)) for x in self.linkdict.keys()], key=lambda x: x[1])[0][0]
                self.update_neighbors(timestep, mode)
            self.results.append({'timeline':[0],
                                'prevalence':[sum(self.state.values())],
                                'ISlinks':[self.countISLinks()],
                                'meta': { 'alpha': self.alpha,
                                        'beta': self.beta,
                                        'R0': self.R0,
                                        'stepsize': self.stepsize,
                                        'infected': [user for user,state in self.state.items() if state],
                                        'duration': self.duration,
                                        'mode': mode
                                }})

            [self.calculate_recovery_process(node) for node,state in self.state.items() if state]


            if mode == 'daynight':
                epsilon = self.starting_epsilon
            while  self.duration > cur_time:
                if epsilon == 1:
                    if mode == 'daynight':
                        self.update_neighbors(cur_time, mode)
                    else:
                        timestep = sorted([(x,abs(cur_time % self.duration_data - x)) for x in self.linkdict.keys()], key=lambda x: x[1])[0][0]
                        self.update_neighbors(timestep, mode)
                self.calculate_infection_process()
                self.update_lambda()
                capLambda= sum(self.s_lambda.values())
                #if sum(self.state.values()) < 285 and cur_time > 1.496 : import pdb; pdb.set_trace()
                #import pdb; pdb.set_trace()
                if epsilon * capLambda * self.stepsize <= tau:
                    #print "case 1"
                    tau = tau - epsilon * capLambda * self.stepsize
                    cur_time += epsilon * self.stepsize
                    epsilon = 1

                else: # epsilon * capLambda * self.stepsize > tau
                    # distribution with probabilities according to rates of the processes (lambda_m / cap_lambda)
                    #print "case 2"

                    prob_dist=rv_discrete(name='cookies',
                                            values=(self.user_idx,
                                                [self.s_lambda[self.users[x]]/capLambda for x in self.user_idx]))
                    node_for_event=self.users[prob_dist.rvs()]
                    cur_time += tau / capLambda

                    self.update_state(node_for_event,cur_time)
                    #it is important that the state was updated BEFORE the recovery process and result update
                    self.calculate_recovery_process(node_for_event)
                    ISlinks = self.countISLinks()
                    self.update_results(node_for_event, ISlinks, cur_time)
                    epsilon -= tau / (capLambda * self.stepsize)
                    tau = np.random.exponential(1)

            counter +=1
        if mode == "daynight":
            self.stepsize = self.regular_stepsize
        return

    def save(self,path,mkdir=True):
        """ writes the simulation results in a pickle file; if mkdir is set on true creates the path if missing; saved results are saved as list of lists"""
        import pickle
        import os

        if mkdir:
            import os
            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))

        pickle.dump(self.results,open(path,'w'))
