def read_file(path_to_file):
    #read file
    with open(path_to_file, 'r') as f:
        link_list = [line.rstrip('\r\n').split(' ') for line in f]
    # remove selflinks
    links=[[x[0],x[1],int(x[2]),int(x[3])] for x in link_list[1:] if not x[0]==x[1]]

    return links


def shuffle_keepdegree_ext(path_to_file, iterations, as_graph = True, return_original = False, write_to_file = False):
    '''
    take data file and shuffle the links while keeping the degree from the original data
    '''

    from random import choice
    import numpy as np
    import sys



    #read file
    with open(path_to_file, 'r') as f:
        link_list = [line.rstrip('\r\n').split(' ') for line in f]
    # remove selflinks
    links=[[x[0],x[1],int(x[2]),int(x[3])] for x in link_list[1:] if not x[0]==x[1]]
    lilen=len(links)
    if not all([x[3]-links[0][3]==0 for x in links]):
        print "not all links have the same duration, unable to shuffle"
        return

    # oridict = {}
    N = len(set([y for x in links for y in [x[0],x[1]]]))
    timeframes=sorted(list(set([x[2] for x in links])))
    timestep = links[0][3]
    timespan = timeframes[-1] + timestep - timeframes[0]
    maxstep=len(timeframes)
    megalist = []

    clustering = {'original': [], 'shuffled': [], 'random': []}

    if return_original: original = []
    for counter,timestamp in enumerate(timeframes):

        sys.stdout.write("\r Step "+str(counter)+' of '+str(maxstep))
        sys.stdout.flush()
        links_for_ts = []
        neighbors_for_ts = {}
        #get all the links for timestep while removing them from the list
        to_delete=[]
        for idx,link in enumerate(links):
            if link[2]==timestamp:
                links_for_ts.append(links[idx])
                to_delete.append(idx)
        for idx in sorted(to_delete,reverse=True): del links[idx]

        #check for duplicate links, remove if present
        first=[]
        link_set=[frozenset([link[0],link[1]]) for link in links_for_ts]
        if len(set(link_set))<len(links_for_ts):
            duplicates=sorted([idx for idx, x in enumerate(link_set) if link_set.count(x) > 1],reversed=True)
            for dup in duplicates:
                if not links_for_ts[dup] in first:
                    first.append(links_for_ts[dup])
                else:
                    del links_for_ts[dup]
        # oridict[timestamp] = [(x[0],x[1]) for x in links_for_ts]
        if return_original: original.extend([(x[0],x[1]) for x in links_for_ts])

        #run shuffling algorithm
        shuffled_list = randomize_by_edge_swaps(links_for_ts, iterations)
        megalist.extend(shuffled_list)
        clustering['shuffled'].append(calculate_clustering_coef(shuffled_list))
        clustering['original'].append(calculate_clustering_coef([(x[0],x[1]) for x in links_for_ts]))
        clustering['random'].append(len(links_for_ts)/(N/2*(N-1)))

        counter += 1

    if write_to_file:
        with open(write_to_file, 'w') as f:
            f.write('# user1 user2 timestamp duration\n')
            for link in megalist:
                f.write(link[0]+' '
                        +link[1]+' '
                        +str(link[2])+' '
                        +str(link[3]))
                f.write('\n')

    # return megadict, oridict
    if as_graph:
        if return_original:
            return list_to_graph(megalist, timespan, timestep, nodes_as_int = True, directed = False), list_to_graph(original, timespan, timestep, nodes_as_int = True, directed = False), clustering
        else: return list_to_graph(megalist, timespan, timestep, nodes_as_int = True, directed = False), clustering

    if return_original: return megalist, original, timeframes
    else: return megalist







#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Copyright (c) 2011-2012 Christopher D. Lasher
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



def randomize_by_edge_swaps(list_of_edges, num_iterations, avoid_repetitions = True):
    """Randomizes the graph by swapping edges in such a way that
    preserves the degree distribution of the original graph.

    The underlying idea stems from the following. Say we have this
    original formation of edges in the original graph:

        head1   head2
          |       |
          |       |
          |       |
        tail1   tail2

    Then we wish to swap the edges between these four nodes as one
    of the two following possibilities:

        head1   head2       head1---head2
              X
        tail1   tail2       tail1---tail2

    We approach both by following through the first of the two
    possibilities, but before committing to the edge creation, give
    a chance that we flip the nodes head1' and 'tail1'.

    :Parameters:
    - 'num_iterations': the number of iterations for edge swapping
      to perform; this value will be multiplied by the number of
      edges in the graph to get the total number of iterations
      """

    import itertools
    import random
    import sys

    from copy import deepcopy
    if type(list_of_edges[0])==list : edge_list = [(x[0],x[1]) for x in list_of_edges]
    else: edge_list = deepcopy(list_of_edges)
    num_edges = len(edge_list)
    total_iterations = num_edges * num_iterations
    ptone_pc = min(max(1,total_iterations / 1000),10000)

    for i in xrange(total_iterations):
        if i % ptone_pc == 0:
            perc = float(i)/total_iterations
            sys.stdout.write("\r"+str(perc)+" % done")
            sys.stdout.flush()
        rand_index1 = int(round(random.random() * (num_edges - 1)))
        rand_index2 = int(round(random.random() * (num_edges - 1)))
        original_edge1 = edge_list[rand_index1]
        original_edge2 = edge_list[rand_index2]
        head1, tail1 = original_edge1[0], original_edge1[1]
        head2, tail2 = original_edge2[0], original_edge2[1]

        # Flip a coin to see if we should swap head1 and tail1 for
        # the connections
        if random.random() >= 0.5:
            head1, tail1 = tail1, head1

        # The plan now is to pair head1 with tail2, and head2 with
        # tail1
        #
        # To avoid self-loops in the graph, we have to check that,
        # by pairing head1 with tail2 (respectively, head2 with
        # tail1) that head1 and tail2 are not actually the same
        # node. For example, suppose we have the edges (a, b) and
        # (b, c) to swap.
        #
        #   b
        #  / \
        # a   c
        #
        # We would have new edges (a, c) and (b, b) if we didn't do
        # this check.

        if head1 == tail2 or head2 == tail1:
            continue

        # Trying to avoid multiple edges between same pair of nodes;
        # for example, suppose we had the following
        #
        # a   c
        # |*  |           | original edge pair being looked at
        # | * |
        # |  *|           * existing edge, not in the pair
        # b   d
        #
        # Then we might accidentally create yet another (a, d) edge.
        # Note that this also solves the case of the following,
        # missed by the first check, for edges (a, b) and (a, c)
        #
        #   a
        #  / \
        # b   c
        #
        # These edges already exist.
        if avoid_repetitions and any([x in edge_list for x in [ (head1, tail2), (tail2, head1), (head2, tail1), (tail1,head2) ] ]):
            continue

        # Suceeded checks, perform the swap
        # update the entries at the indices randomly selected
        edge_list[rand_index1] = (head1, tail2)
        edge_list[rand_index2] = (head2, tail1)

    assert len(edge_list) == num_edges
    return edge_list


def list_to_graph(edge_list, timespan, timestep, nodes_as_int = False, users = None, exportIDdict = False):
    import networkx as nx

    G = nx.Graph()

    if nodes_as_int:

        if type(nodes_as_int)==dict: ID_to_int = nodes_as_int
        else: ID_to_int = {x : idx for idx,x in enumerate(set([y for x in edge_list for y in [x[0],x[1]]]))}
        for node_name in ID_to_int.keys():
            G.add_node(ID_to_int[node_name], name=node_name)
        for edge in edge_list:
            if G.has_edge(ID_to_int[edge[0]],ID_to_int[edge[1]]): G.edge[ID_to_int[edge[0]]][ID_to_int[edge[1]]]['weight'] += 1
            elif G.has_edge(ID_to_int[edge[1]],ID_to_int[edge[0]]): G.edge[ID_to_int[edge[1]]][ID_to_int[edge[0]]]['weight'] += 1
            else: G.add_edge(ID_to_int[edge[0]],ID_to_int[edge[1]], weight = 1)
    else:
        if users: [G.add_node(x) for x in users]
        for edge in edge_list:
            if G.has_edge(edge[0],edge[1]): G.edge[edge[0]][edge[1]]['weight'] += 1
            elif G.has_edge(edge[1],edge[0]): G.edge[edge[1]][edge[0]]['weight'] += 1
            else: G.add_edge(edge[0],edge[1], weight = 1)

    timespan /= float(timestep)
    G = nx.DiGraph(G)
    for edge in G.edges(data=True):
        G.add_edge(edge[0],edge[1], weight = edge[2]["weight"]/timespan)
    if exportIDdict: return G, ID_to_int
    else: return G

def list_to_daynight(edge_list,interval, nodes_as_int = False, exportIDdict=False):

    from datetime import datetime as dt

    timestep = edge_list[0][3]
    day = []
    night = []
    #sorting
    for edge in edge_list:
        linktime=dt.fromtimestamp(float(edge[2]))
        if linktime.hour<interval[0] or linktime.hour>=interval[1]:
            night.append(edge)
        else:
            day.append(edge)

    #boundries
    day_duration = [day[0][2]]
    for idx,x in enumerate(day[:-1]):
        if day[idx+1][2] - x[2] > timestep:
            day_duration.append(x[2]+timestep)
            day_duration.append(day[idx+1][2])
    day_duration.append(day[-1][2]+timestep)
    day_duration = sum([day_duration[x+1] - day_duration[x] for x in range(0,len(day_duration),2)])

    night_duration = [night[0][2]]
    for idx,x in enumerate(night[:-1]):
        if night[idx+1][2] - x[2] > timestep:
            night_duration.append(x[2]+timestep)
            night_duration.append(night[idx+1][2])
    night_duration.append(night[-1][2]+timestep)
    night_duration = sum([night_duration[x+1] - night_duration[x] for x in range(0,len(night_duration),2)])

    if nodes_as_int:
        IDdct = {x : idx for idx,x in enumerate( set([y for x in edge_list for y in [x[0],x[1]]]) ) }
        G_day = list_to_graph(day, day_duration, timestep, nodes_as_int = IDdct, exportIDdict = False)
        G_night = list_to_graph(night, night_duration, timestep, nodes_as_int = IDdct, exportIDdict = False)
    else:
        user = set([y for x in edge_list for y in [x[0],x[1]]])
        G_day = list_to_graph(day, day_duration, timestep, users = user, exportIDdict = False)
        G_night = list_to_graph(night, night_duration, timestep, users = user, exportIDdict = False)
    if exportIDdict: return [G_day,G_night,IDdct]
    else: return [G_day,G_night]


def calculate_clustering_coef(edge_list):
    from networkx import Graph, clustering
    from scipy import mean
    G = Graph()
    G.add_edges_from(edge_list)
    return mean(clustering(G).values())
