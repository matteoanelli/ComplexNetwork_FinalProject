import numpy as np
import networkx as nx
import scipy
from matplotlib import pyplot as pl
import si_animator as an
import pandas as pd
import random
def SI(seed, p, data):
    infection_times = {}
    infection_list = []
    infection_times[seed] = 1229231100
    for index, event in data.iterrows():
        if event['Source'] in infection_times and infection_times[event['Source']] <= event['StartTime']:
            if event['Destination'] not in infection_times:
                if p >= random.randrange(0,1):
                    infection_times[event['Destination']] = event['EndTime']
            else:
                if event['EndTime'] < infection_times[event['Destination']]:
                    infection_times[event['Destination']] = event['EndTime']
    infection_list = list(infection_times.values())

    return infection_times, sorted(infection_list)

def main():
    # event_data = np.genfromtxt('events_US_air_traffic_GMT.txt', names=True, dtype=int)
    # print(event_data)
    # sorted_data = event_data.sort(order=['StartTime'])
    # print(sorted_data)

    data = pd.read_csv('events_US_air_traffic_GMT.txt', header=[0,1,2,3,4], sep=' ')
    data.columns = ['Source', 'Destination', 'StartTime', 'EndTime', 'Duration']
    sorted_data = data.sort_values(by=['StartTime'])
    infection_times, infection_list = SI(0,1,sorted_data)
    print(infection_times[41])

    # task 2

    # prob = [0.01, 0.05, 0.1, 0.5, 1.0]
    # infection_results = {}
    # for p in prob:
    #     infection_list_temp = []
    #     for _ in range(10):
    #         infection_times, infection_list = SI(0,p,sorted_data)
    #         infection_list_temp.append()










if __name__ == '__main__':
    main()