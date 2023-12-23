'''
Title: Multi-scale Multipopulation Epidemic Simulation with the SIRS model.
Description: The simulation below uses a graphing library 'networkx' to simulate
              the evolution of an epidemic over a global population that is comprised
              of 4 local subpopulations. At the local level, the model assumes random
              mixing and applies the SIRS epidemic model. At the global level, it allows
              individuals to move from one population to another, simulating the notion of 
              shortcuts in a small-world network.
Author: Eugene Halim
'''
# Imports
import networkx as nx
import matplotlib.pyplot as plt
import random

# Ask users for input
gamma = input("Contact probability % (for 0.40 = 40% type '40'): ")
while (not gamma.isdigit() or int(gamma) > 100 or int(gamma) < 0):
  print("Please only input integer values between 0 and 100!")
  gamma = input("Contact probability % (for 0.40 = 40% type '40'): ")

beta = input("Infection probability % (for 0.40 = 40% type '40'): ")
while (not beta.isdigit() or int(beta) > 100 or int(beta) < 0):
  print("Please only input integer values between 0 and 100!")
  beta = input("Infection probability % (for 0.40 = 40% type '40'): ")
  
recovery = input("Recovery rate % (for 0.40 = 40% type '40'): ")
while (not recovery.isdigit() or int(recovery) > 100 or int(recovery) < 0):
  print("Please only input integer values between 0 and 100!")
  recovery = input("Recovery rate % (for 0.40 = 40% type '40'): ")

move_probability = input("Probability (%) that a node moves to new context (for 0.40 = 40% type '40'): ")
while (not move_probability.isdigit() or int(move_probability) > 100 or int(move_probability) < 0):
  print("Please only input integer values between 0 and 100!")
  move_probability = input("Probability (%) that a node moves to new context (for 0.40 = 40% type '40'): ")

n = input("Population (integer) of each local context. Pick number between 0-50: ")
while (not n.isdigit() or int(n) > 50 or int(n) < 0):
  print("Please only input integer values between 0 and 50!")
  n = input("Population (integer) of each local context. Pick number between 0-50: ")

num_iter = input("How many iterations (integer):  ")
while (not num_iter.isdigit() or int(num_iter) < 0):
  print("Please only input integer values greater 0!")
  num_iter = input("How many iterations (integer):  ")

# Convert to integers
gamma = float("0." + gamma)
beta = float("0." + beta)
recovery = float("0." + recovery)
move_probability = float("0." + move_probability)
n = int(n)
num_iter = int(num_iter)
r_to_s = recovery

# Global Variables
graphs = [] # Array of dictionaries, one for each of the subpopulation graphs.
global_susceptible = []
global_infected = []
global_recovered = []
timesteps = []

'''
FUNCTION DEFINITIONS
'''
def generateErdosRenyi(susceptible, infected, recovered, p):
  '''
  Function that generates a random graph. It first adds a total of N = susceptible + infected + recovered
  nodes into the graph (with appropiate labels), then it randomly adds edegs between all pairs of nodes 
  with probability p. 
  Inputs:
  - Susceptible: Integer, number of susceptible individuals in subpopulatoin.
  - Infected: Integer, number of infected individuals in subpopulatoin.
  - Recovered: Integer, number of recovered individuals in subpopulatoin.
  - p: Float between 0 and 1, contact probability between individuals in the same graph. 

  Outputs:
  - G: The networkx graph object.
  - Colormap: Colors to assign to each individual node in graph. Creates visual distinction between different
              node states.
  '''

  # Find total populatoin in graph
  N = susceptible + infected + recovered
  # Create an empty graph
  G = nx.Graph()

  color_map = []

  # Iterate through each node, adding labels and colors 
  for i in range(N):
    if i < susceptible:
      G.add_node(i, label="susceptible")
      color_map.append('orange')

    elif i < susceptible + infected: 
      G.add_node(i, label="infected")
      color_map.append('red')      
    else:
      G.add_node(i, label="recovered")
      color_map.append('green')

  # iterate between all pairs of nodes (i, j)
  for i in range(0, N):
    for j in range(i+1, N):
      # pick a random number and if it is less than p add an edge
      if random.random() < p:
        G.add_edge(i, j)

  # Return the graph with edges and color map
  return G, color_map

def findContextToMoveTo(group_num):
  '''
  Finds which subpopulation an individual from a given group will move to. The sampling 
  function is a simplified version of the actual multiscale metapopulatin model - groups 
  that are a distance x from the current group are chosen with probability 1/2^x. The 
  mapping below was written or a global population with only 4 subpopulaton, where the 
  binary tree is of height 2. The following is the structure of the binary tree:

                    /\
                   /  \
                  /\  /\
                 1 2 3 4
  Inputs:
  - group_num: The subgroup that the individual is currently in.

  Outputs:
  - Integer: The subgroup that the given individual is moving to. 
  '''

  # Generate a random number
  sampler = random.random()

  # With probability 50%, return the group that is distance 1 away. 
  if (sampler <= 0.5):
    if (group_num == 1):
      return 2
    elif (group_num == 2):
      return 1
    elif (group_num == 3):
      return 4
    elif (group_num == 4):
      return 3
    
  # With probability 25%, return a group that is distance 2 away.
  elif(sampler > 0.5 and sampler <=0.75):
    if (group_num == 1):
      return 3
    elif (group_num == 2):
      return 3
    elif (group_num == 3):
      return 1
    elif (group_num == 4):
      return 1
    
  # With probability 25%, return the second group tat is distance 2 away.
  else:
    if (group_num == 1):
      return 4
    elif (group_num == 2):
      return 4
    elif (group_num == 3):
      return 2
    elif (group_num == 4):
      return 2
   

def performGraphIteration(g, num_graph):
  '''
  Runs the simulation for a single timestep for the given graph. Iterates through each node in the 
  graph and decides the next state of the current node based on its current state and the probability
  of state transitions that the user inputs. It updates both the local and global population statistics.
  This function does not redraw the graph, but provides the necessary information for the "generateErdosRenyi"
  function to generate the corresponding graph.

  Inputs:
  - g: the NetworkX graph object
  - num_graph: The graph's position in the global list of graph dictionaries. 

  Outputs: None. Performs updates in-place.
  '''

  # Get the statistics from the last timestep from the current graph's disctionary.
  local_susceptible = graphs[num_graph]["susceptible"]
  local_infected = graphs[num_graph]["infected"]
  local_recovered = graphs[num_graph]["recovered"]

  # Nodes to remove in case a node moves from the current subgroup to another.
  nodes_to_remove = []

  # iterate through each node
  for node in g:
    # Check if the node will move to a new subpopulation. 
    if (random.random() < move_probability):
      new_context = findContextToMoveTo(num_graph + 1)

      # Add the current node with the current label to new graph and update population count.
      graphs[new_context - 1]["graph"].add_node(node, label=g.nodes[node]['label'])
      graphs[new_context - 1][g.nodes[node]['label']] = graphs[new_context - 1][g.nodes[node]['label']] + 1

      # Save node information to remove later.
      nodes_to_remove.append(node)

    # If infected, check if node recovers within current timestep.
    elif (g.nodes[node]['label'] == "infected" and random.random() < recovery): 
      # Update node label and graph population statistics
      g.nodes[node]['label'] = "recovered"
      local_infected = local_infected - 1
      local_recovered = local_recovered + 1

    # If recovered, check whether becomes susceptible again.
    elif (g.nodes[node]['label'] == "recovered" and random.random() < r_to_s):
      # Make current node susceptible and update population statistics.
      g.nodes[node]['label'] = "susceptible"
      local_recovered = local_recovered - 1
      local_susceptible = local_susceptible + 1

    # If susceptible, check if infected by neighbors based on probability given
    elif (g.nodes[node]['label'] == "susceptible"):

      # Iterate through each of the neighboring nodes.
      newly_infected = []
      for neighbor in g[node]:
        # If neighbor is infected, sample probability that current node gets infected.
        if (g.nodes[neighbor]['label'] == "infected" and random.random() < beta):
          # Change current node's state if probability falls in given threshold 
          newly_infected.append(node)
          local_susceptible = local_susceptible - 1
          local_infected = local_infected + 1

          # break since every node can only get infected once. Prevents double counting.
          break
    
      # Only update states after timestep since nodes only become infected after 1 timestep.
      for x in newly_infected:
        g.nodes[x]['label'] = "infected"

  # Only remove node after timestep to avoid errors.
  for node in nodes_to_remove:
    # First update population information
    if (g.nodes[node]['label'] == "susceptible"):
      local_susceptible = local_susceptible - 1
    elif (g.nodes[node]['label'] == "infected"):
      local_infected = local_infected - 1
    elif (g.nodes[node]['label'] == "recovered"):
      local_recovered = local_recovered - 1
    
    # Then remove node.
    g.remove_node(node)

  # update graph dictionary:
  graphs[num_graph]["susceptible"] = local_susceptible
  graphs[num_graph]["infected"] = local_infected
  graphs[num_graph]["recovered"] = local_recovered
      
'''
MAIN CODE
'''
# Main while loop
t = 0
while t <= num_iter:
  if (t == 0):
    # Generate initial conditions for each graph (1 infected each. The rest are susceptible)
    for i in range(4):
      (g, cmap) = generateErdosRenyi(n-1, 1, 0, gamma)
      graph = {
        "graph": g, 
        "cmap": cmap, 
        "susceptible": n-1,
        "infected": 1,
        "recovered" : 0
        }
      graphs.append(graph)

    # Update global states:
    global_infected.append(4)
    global_recovered.append(0)
    global_susceptible.append(4 * (n-1))
    timesteps.append(0)

  else: # For every other timestep:
    # Keep track of statistics for current timestep
    curr_global_susceptible = 0
    curr_global_infected = 0
    curr_global_recovered = 0

    for n in range(len(graphs)):
      # Perform 1 timestep of epidemic spread
      performGraphIteration(graphs[n]["graph"], n)

      # Extract total S,I,R populations after updating local graphs
      curr_global_susceptible = curr_global_susceptible + graphs[n]["susceptible"]
      curr_global_infected = curr_global_infected + graphs[n]["infected"]
      curr_global_recovered = curr_global_recovered + graphs[n]["recovered"]

    # Redraw graphs after all state updates complete to prevent undrawable graphs:
    for n in range(len(graphs)):
      # Regenerate graph edges randomly (with same states) to simulate random mixing
      (g, cmap) = generateErdosRenyi(graphs[n]["susceptible"], graphs[n]["infected"], graphs[n]["recovered"], gamma)

      # Save new graphs to graph dictionary
      graphs[n]["graph"] = g
      graphs[n]["cmap"] = cmap

    # Update global variables
    global_susceptible.append(curr_global_susceptible)
    global_infected.append(curr_global_infected)
    global_recovered.append(curr_global_recovered)
    timesteps.append(t)

  # Draw graphs:
  plt.figure(figsize=(10,8))
  for n in range(len(graphs)):
    population = graphs[n]["susceptible"] + graphs[n]["infected"] + graphs[n]["recovered"] 
    plt.subplot(int("32" + str(n+1)))
    plt.title("Graph" + str(n+1) + ", iter=" + str(t) + " , n = " + str(population), fontsize = 10)
    plt.xlabel("S = " + str(graphs[n]["susceptible"]) + ", I = " + str(graphs[n]["infected"]) + ", R = " + str(graphs[n]["recovered"]), fontsize=9)
    nx.draw_networkx(graphs[n]["graph"], node_color=graphs[n]["cmap"])

  # Plot global statistics:
  plt.subplot(int("32" + str(5)))
  plt.plot(timesteps, global_susceptible, "-b", label="susceptible")
  plt.plot(timesteps, global_infected, "-r", label="infected")
  plt.plot(timesteps, global_recovered, "-g", label="recovered")
  plt.title("Global Statistics")
  plt.xlabel("timestep")
  plt.ylabel("# people")
  plt.legend()

  plt.show()

  # Increment timestep
  t = t + 1