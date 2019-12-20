import numpy as np
import networkx as nx
import scipy
from scipy.stats import spearmanr
from matplotlib import pyplot as plt
from collections import defaultdict, Counter
import si_animator as an
import pandas as pd
from itertools import groupby
import random
from si_animator import visualize_si, plot_network_usa
from tqdm import tqdm

def SI(seed, p,data,time_inervals=None, immune_nodes = []):
    infection_times = {}

    # Set starting node and time
    infection_times[seed] = 1229231100

    for event in data:
        if event['Source'] in infection_times and infection_times[event['Source']] <= event['StartTime'] and\
            event['Destination'] not in immune_nodes:
            if event['Destination'] not in infection_times:
                if random.random() <= p:
                    infection_times[event['Destination']] = event['EndTime']
            else:
                if event['EndTime'] < infection_times[event['Destination']]:
                    infection_times[event["Destination"]] = event['EndTime']

    infection_list = list(infection_times.values())
    infection_list.sort()

    if time_inervals is None:
        return infection_times, infection_list
    else:
        return infection_times, list(np.digitize(infection_list, time_inervals))
def average_prevelance(infection_list, bins, number_nodes):
    average_list = []
    total = 0
    for bin in range(bins):
        if bin+1 in infection_list:
            count = infection_list.count(bin+1)
            total += count
            average_list.append(total /number_nodes)
        else:
            average_list.append(total/number_nodes)
    return average_list

def infection_multiple(probs, times, infected_source, sorted_data, bins, n_bins, number_nodes):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for p in probs:
        list_prob_times = []
        for _ in range(times):
            _, infection_list = SI(infected_source, p, sorted_data, bins)
            list_times = average_prevelance(infection_list, n_bins, number_nodes)
            list_prob_times.append(list_times)
        plt.plot(bins-1229231100, np.mean(list_prob_times, 0),
                 '-',
                 linewidth=1.0,
                 label='p={:.2f}'.format(p)
                 )

    ax.set_title('Average prevalence as a function of time different probabilities')
    ax.set_ylabel(r'Average prevalence')
    ax.set_xlabel(r'Time')
    ax.legend()

    return fig

def infection_multiple_seed(p, times, seeds, sorted_data, bins, n_bins, number_nodes):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    aereport = {0: 'ABE',
                4: 'ATL',
                41: 'ACN',
                100: 'HSV',
                200: 'DBQ'}

    for seed in seeds:
        list_prob_times = []
        for _ in range(times):
            _, infection_list = SI(seed, p, sorted_data, bins)
            list_times = average_prevelance(infection_list, n_bins, number_nodes)
            list_prob_times.append(list_times)
        plt.plot(bins-1229231100, np.mean(list_prob_times, 0),
                 '-',
                 linewidth=1.0,
                 label=aereport[seed]
                 )

    ax.set_title('Average prevalence as a function of seed different seeds')
    ax.set_ylabel(r'Average prevalence')
    ax.set_xlabel(r'Time')
    ax.legend('lower right')

    return fig

# def order_dict(values):
#     temp = [val for val in list(values) if val[1]]
#     temp.sort(key=lambda x: x[1])
#     return temp

def median_infection(network: nx.Graph, sorted_data, p=.5, times=20):
    dict_infection = defaultdict(list)
    for _ in range(times):
        random_seed = int(random.choice(list(network.nodes())))
        infection_times, _ = SI(random_seed, p, sorted_data)
        for key, value in infection_times.items():
            if key in dict_infection:
                dict_infection[key].append(value)
            else:
                dict_infection[key].append(value)

    median_dict = {}
    for node in dict_infection:
        median_dict[node] = np.median(dict_infection[node])

    return median_dict

def similarity_measures(network: nx.Graph, sorted_data, p=.5, times=50):
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    median_infection_dict = median_infection(network,sorted_data)

    similarity_measure_dict = {'UCC' : nx.clustering(network),
                               'Degree': nx.degree(network),
                               'Strenght': nx.degree(network,weight='weight'),
                               'UBC': nx.betweenness_centrality(network)}

    for key, value in similarity_measure_dict.items():
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)

        # Get for each node the measure value
        x = [value[node] for node in list(network.nodes)]
        # get the median for each node
        y = [median_infection_dict[int(node)] for node in list(network.nodes)]

        ax.set_title(r'Median infection time as a function of {}'.format(key))
        ax.set_ylabel(r'{}'.format(key))
        ax.set_xlabel(r'Median infection time')
        if key != 'UCC' and key != 'UBC':
            ax.set_xscale('log')

        plt.scatter(x,y)
        fig.savefig('./task4_{}.pdf'.format(key))
        fig.show()
        fig.clear()

        spearman_coeff = spearmanr(x, y).correlation
        print('Spearman rank-correlation coefficient of {}: {:.4f}'.format(key, spearman_coeff))

def random_neighbors(network, number_of_nodes):
    random_neighbors = []
    while len(random_neighbors) < number_of_nodes:
        neighbors = nx.neighbors(network, random.choice(list(network.nodes)))
        random_neighbor = random.choice(list(neighbors))
        if random_neighbor not in random_neighbors:
            random_neighbors.append(random_neighbor)
    return random_neighbors

def immune_nodes(network :nx.Graph, number_of_nodes):
    strategies = ["random neighbor","random node","clustering coefficient","degree","strength","betweenness centrality"]

    immune_nodes = {
        'Social Net.': random_neighbors(network, number_of_nodes),
        'Random Node': random.sample([int(n) for n in list(network.nodes)],number_of_nodes),
        'UCC': [int(i[0]) for i in Counter(nx.clustering(network)).most_common(number_of_nodes)],
        'Degree': [int(i[0]) for i in Counter(dict(nx.degree(network))).most_common(number_of_nodes)],
        'Strenght': [ int(i[0][0])for i in Counter((nx.degree(network,weight='weight'))).most_common(number_of_nodes)],
        'UBC': [int(i[0]) for i in Counter(nx.betweenness_centrality(network)).most_common(number_of_nodes)]}
    return immune_nodes

def shutting_down_airports(network :nx.Graph, sorted_data, immune_nodes,bins, n_bins,number_of_nodes, p=.5, times=20):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)

    set_immune_node = []
    for key in immune_nodes:
        set_immune_node += immune_nodes[key]


    seeds = []

    while len(seeds) < 20:
        node = int(random.choice(list(network.nodes)))
        if node not in set_immune_node:
            seeds.append(node)

    for strategy in immune_nodes.keys():
        prevalence = []
        for seed in  tqdm(seeds):
            _, infection_list = SI(seed, .5, sorted_data, bins,immune_nodes[strategy])
            list_times = average_prevelance(infection_list, n_bins, number_of_nodes)
            prevalence.append(list_times)

        plt.plot(bins - 1229231100, np.mean(prevalence, 0),
                 '-',
                 linewidth=1.0,
                 label='{}'.format(strategy)
                 )

    ax.set_title(r'Average prevalence as a function of time with different immunity strategy')
    ax.set_ylabel(r'Average prevalence')
    ax.set_xlabel(r'Time')
    ax.legend(loc=0)

    return fig

def SI_links(sorted_data, seed, p=.5):
    infection_times = {}
    # Set starting node and time
    infection_times[seed] = 1229231100
    links ={}

    for event in sorted_data:
        if event['Source'] in infection_times and infection_times[event['Source']] <= event['StartTime']:
            if event['Destination'] not in infection_times:
                if random.random() <= p:
                    # Before replace the time than create the infected link
                    infection_times[event['Destination']] = event['EndTime']
                    links[(event['Source']), event['Destination']] = 1
            else:
                if event['EndTime'] < infection_times[event['Destination']]:
                    infection_times[event["Destination"]] = event['EndTime']
                    # increase the number of the specific infected link
                    links.pop(([ link for link in links.keys() if link[1] == event["Destination"]][0]))
                    links[(event["Source"], event["Destination"])] = 1

    return links


def compute_infection_links(network, event_data, p=.5, times=20):
    nodes = list(network.nodes)

    # random seed to start
    seeds = random.sample([int(n) for n in nodes], times)

    links_dict = {}
    list_dict = []
    for seed in seeds:
         list_dict.append(SI_links(event_data, seed))

    for d in list_dict:
        for key, value in d.items():
            if key not in links_dict.keys():
                links_dict[key] = 1
            else:
                links_dict[key] += 1

     # order airport key a<b

    links_dict_ordered = {}
    for key, value in links_dict.items():
        if (key[1], key[0]) in links_dict:
            if (str(min(key[0], key[1])), str(max(key[0], key[1]))) in links_dict_ordered:
                links_dict_ordered[str(min(key[0], key[1])), str(max(key[0], key[1]))] += value
            else:
                links_dict_ordered[str(min(key[0], key[1])), str(max(key[0], key[1]))] = value
        else:
            links_dict_ordered[str(min(key[0], key[1])), str(max(key[0], key[1]))] = value
    return links_dict_ordered

def scatter_links(network, weights):
    edges = list(network.edges())

    statistics = {
        'Weight': nx.get_edge_attributes(network, 'weight'),
        'Link Betweenness Centrality': nx.edge_betweenness_centrality(network)
    }
    # fill all the weights also for edges that are missing

    weights_full = []
    for edge in edges:
        if edge not in weights:
            weights_full.append(0)
        else:
            weights_full.append(weights[edge]/20)

    for key, value in statistics.items():
        x = [value[edge] for edge in edges]
        y = weights_full

        plt.xlabel('{}'.format(key))
        plt.ylabel('Fraction of times fij')
        plt.scatter(x, y)
        plt.title('Transmission ratio as function of {}'.format(key))
        plt.savefig('./{}.pdf'.format(key.replace(" ", "_")))
        plt.show()
        plt.close()

        spearman_coeff = spearmanr(x, y).correlation
        print('Spearman rank-correlation coefficient of {}: {:.4f}'.format(key, spearman_coeff))




def main():

    ndtype = [("Source", int), ("Destination", int), ("StartTime", int), ("EndTime", int)]
    event_data = np.genfromtxt('events_US_air_traffic_GMT.txt', names=True, dtype=ndtype)
    event_data.sort(order=["StartTime"])


    # reading input event data and sort them in increasing departure time

    # data = pd.read_csv('events_US_air_traffic_GMT.txt', header=[0,1,2,3,4], sep=' ')
    # data.columns = ['Source', 'Destination', 'StartTime', 'EndTime', 'Duration']
    # sorted_data = data.sort_values(by=['StartTime'])

    # Reading network
    with open('./aggregated_US_air_traffic_network_undir.edg', 'rb') as edge_list:
        network = nx.read_weighted_edgelist(edge_list)

    data = np.genfromtxt("./US_airport_id_info.csv", delimiter=',', dtype=None, names=True, encoding=None)
    xycoords = {}
    for row in data:
        xycoords[str(row['id'])] = (row['xcoordviz'], row['ycoordviz'])

    nodes = network.nodes
    number_nodes = len(nodes)

    # fix number of bins
    n_bins = 75
    bins = np.linspace(1229231100, 1230128400, 75).astype(int)

    # ----------------------------------------------------- Task 1 -----------------------------------------------------#
    print('Task1........')
    infection_times, infection_list = SI(0,1,event_data)
    print(infection_times[41])
    print('Task1 end')
    #----------------------------------------------------- Task 2 -----------------------------------------------------#
    print('Task2........')
    prob = [0.01, 0.05, 0.1, 0.5, 1.0]

    fig = infection_multiple(prob, 10, 0, event_data, bins, n_bins, number_nodes)
    fig.savefig('./task2_average_prevalence.pdf')
    fig.show()
    fig.clear()
    print('Task2 end')


    #----------------------------------------------------- Task 3 -----------------------------------------------------#
    print('Task3........')

    seed = [0, 4, 41, 100, 200]

    fig = infection_multiple_seed(.1, 10, seed, event_data, bins, n_bins, number_nodes)
    fig.savefig('./task3_seed_inspection.pdf')
    fig.show()
    fig.clear()
    print('Task3 end')

    #----------------------------------------------------- Task 4 -----------------------------------------------------#
    print('Task4........')

    similarity_measures(network,event_data)
    print('Task4 end')
    #----------------------------------------------------- Task 5 -----------------------------------------------------#
    print('Task5.......')
    number_immune_nodes = 10

    immune_nodes_dict = immune_nodes(network, number_immune_nodes)

    fig = shutting_down_airports(network, event_data, immune_nodes_dict,bins, n_bins,number_nodes)
    fig.savefig('./task5_immunization_analysis.pdf')
    fig.show()
    fig.clear()
    print('Task5 end')

    # ----------------------------------------------------- Task 6 -------------------------------------------- #
    print('Task6........')

    times = 20
    seeds = []
    while len(seeds) < 20:
        node = int(random.choice(list(network.nodes)))
        if node not in seeds:
            seeds.append(node)

    weights = compute_infection_links(network, event_data)
    norm_weights = []
    for l in list(weights.values()):
        norm_weights.append(float(l/20))

    plot_network_usa(network, xycoords, edges=list(weights.keys()), linewidths=norm_weights)
    plt.title("Infected links")
    plt.savefig('./task6_plot_network_usa.pdf')
    plt.show()
    plt.close()

    max_spann_tree = nx.maximum_spanning_tree(network)
    max_spann_edges = list(max_spann_tree.edges)

    plot_network_usa(max_spann_tree, xycoords, edges=max_spann_edges)
    plt.title("Infected links max spanning")
    plt.savefig('./task6_plot_network_usa_max_spanning_tree.pdf')
    plt.show()
    plt.close()

    scatter_links(network, weights)

    print('Task6 end')










if __name__ == '__main__':
    main()