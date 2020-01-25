### Complex Network Final Project
Final project of complex network. (CS-E5740)

## Author 
Matteo Anelli

## Introduction
In the project work, your task is to implement a Susceptible-Infected (SI) disease spreading
model and run it on top of a temporal network from air transport data, containing information
on the departure and arrival times of flights. You will study the dynamics of spreading and
how it depends on where the process starts as well as the infectivity of the disease, and use
static-network centrality measures to understand the roles that specific nodes play.

## Model Specfications 

In the SI model, each node is either Susceptible or Infected. When an Infected node is in
contact with a Susceptible node, the Susceptible node may become infected with some probability
p ∈ [0, 1], reflecting the infectivity of the disease. Infected nodes remain Infected forever.
In our model that mimics the spreading of disease through the air transport network, nodes are
airports and time-stamped connections are flights between them. Initially, only one airport node
(called the seed node) is set to the Infected state, and all other airports are susceptible. Now,
following the SI process, a Susceptible node may become Infected when an airplane originating
from an Infected node arrives to the airport. Note that the source airport needs to be Infected
already at the time of the flight’s departure! If this condition is met, the destination airport
may become Infected with probability p ∈ [0, 1], and if so, it remains in this state for the rest of
the simulation.

## Data description 
The flight data that you will use to perform your simulation are located in the file
__events_US_air_traffic_GMT.txt__, where each row contains the following fields:
- Source [0-278]
- Destination [0-278]
- Start Time (GMT) [seconds after Unix epoch time]
- End Time (GMT)
- Duration [Same as (EndTime-StartTime)]
The aggregated weighted network __aggregated_US_air_traffic_network_undir.edg__ is con-
structed based on the event data, so that the weight of each link corresponds to the number of
flights between the nodes. This network is static and undirected. Additionally, you can find the information about the airports in file __US_airport_id_info.csv__.




