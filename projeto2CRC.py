# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:17:22 2020

@author: Kikoze
"""
#%matplotlib inline
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd     
import random as rd

def corr(total, g, viz, tot, d):
    # Calculate correlation
    for i in range(total):  
        for j in range(total):  
            for p in range(d):
                if (dist1[i][j] == (p+1) and g.nodes[i]['state'] != 0):
                    tot[i][p] += 1
                    if (g.nodes[j]['state'] != 0 and i!=j):
                        viz[i][p] += 1 # count dist1

def states(susceptible, infected, recovered):
    # Evolution in states
    plt.plot(recovered, label = 'Recovered')
    plt.plot(infected, label = 'Infected')
    plt.plot(susceptible, label = 'Susceptible')
    plt.legend()
    plt.xlabel('NRUNS')
    plt.ylabel('N. Nodes')
    plt.show()    
    
# declarar variáveis do modelo SIR
count = 20 
total = 900
k = 6
d = 6
alpha = 0.5
gamma = 1
tempo = 2 # número de iterações até nó poder ficar recuperado
med = [0]*d
medg = [0]*d
medr = [0]*d
for c in range(count):
    
    S = total # Blue
    I = 0 # Red
    R = 0 # Green
    ites = [0]*total
    ite = 0
    susceptible = []
    infected = []
    recovered = []
    
    # Create WS graph
    g = nx.watts_strogatz_graph(total, k, 1, seed=None)
    g_pos = nx.spring_layout(g)
    
    lengths1 = dict(nx.shortest_path_length(g)) 
    dist1 = pd.DataFrame(lengths1) # matriz com as distâncias entre cada ponto
    
    # Simulate SIR Model
    while True:
        if ((I+R)>=(0.5*total)): 
            break
        if(ite != 0 and I == 0):
            break
        if (ite==0):
            # Initialize graph with only one infected
            for i in range(total):
                g.nodes[i]['state'] = 0 
            g.nodes[0]['state'] = 1
            S -= 1
            I += 1
        else:
            for i in range(total):  
                for j in range(total):  
                    if (dist1[i][j] == 1 and g.nodes[j]['state'] == 1 and g.nodes[i]['state'] == 0) and i!=j:
                        aux_alpha = rd.random()
                        if aux_alpha < alpha:
                            g.nodes[i]['state'] = 3 # Susceptible node becomes infected
                            S -= 1
                            I += 1
                            if ((I+R)==(0.5*total)):
                                break
                if ((I+R)==(0.5*total)):
                                break
            for i in range(total): 
                if (g.nodes[i]['state'] == 3):
                    g.nodes[i]['state'] = 1  
                    ites[i] = ite
                
                if (g.nodes[i]['state'] == 1 and ite >= (ites[i]+tempo)):
                    aux_gamma = rd.random()
                    if aux_gamma < gamma:
                        g.nodes[i]['state'] = 2 # Infected node becomes recovered
                        I -= 1
                        R += 1 
                        if ((I+R)==(0.5*total)):
                            break     
            
        # Vector to build the graph temporal evolution
        susceptible.append(S)
        infected.append(I)
        recovered.append(R)
        ite += 1
    
        # Draw graph
        color_map = []
        for node in g:
            if g.nodes[node]['state'] == 0:
                color_map.append('blue')
            elif g.nodes[node]['state'] == 1 or g.nodes[node]['state'] == 3:
                color_map.append('red')
            else: 
                color_map.append('green') 
        
    # Show graph that result from SIR Model simulation
    #nx.draw(g, node_color=color_map, with_labels=True)
    #plt.show()
    
    #print('State S: ',S,'. State I: ',I,'. State R: ',R)
    
    # Show graph temporal evolution plot
    #if(c == 0):
        #states(susceptible, infected, recovered)
    
    # matrix containing number of infected neighbors at distances 
    # d from every infected node from the SIR graph
    gviz = [[0 for x in range(d)] for y in range(total)] 
    # matrix containing number of neighbors at distances 
    # d from every infected node from the SIR
    gtot = [[0 for x in range(d)] for y in range(total)] 
    rv1 = [0]*d #
    rt1 = [0]*d
    r1 = [0]*d
    # Calculate viz and tot vectors from SIR Model graph 
    corr(total, g, gviz, gtot, d)
    
    for i in range(d):
        for j in range(total):
            rv1[i] += gviz[j][i]
            rt1[i] += gtot[j][i]
        r1[i] = rv1[i]/rt1[i]
    
    G = nx.watts_strogatz_graph(total, k, 1, seed=None)
    G_pos = nx.spring_layout(G)
    
    lengths1 = dict(nx.shortest_path_length(G)) 
    
    dist1 = pd.DataFrame(lengths1)
    color_map = []
    for node in G:
        rand = rd.random()
        if rand <= 0.5:
            G.nodes[node]['state'] = 0     
            color_map.append('blue')
        else:
            G.nodes[node]['state'] = 1  
            color_map.append('red')
    
    # Show randomly generated graph
    #nx.draw(G, node_color=color_map, with_labels=True)
    #plt.show()
    
    # matrix containing number of infected neighbors at distance 
    # d from every infected node from the random graph
    rviz = [[0 for x in range(d)] for y in range(total)] 
    # matrix containing number of neighbors at distances 
    # d from every infected node from the random graph
    rtot = [[0 for x in range(d)] for y in range(total)] 
    rv2 = [0]*d
    rt2 = [0]*d
    r2 = [0]*d
    # Calculate viz and tot vectors from random graph 
    corr(total, G, rviz, rtot, d)
    
    for i in range(d):
        for j in range(total):
            rv2[i] += rviz[j][i]
            rt2[i] += rtot[j][i]
        r2[i] = rv2[i]/rt2[i]
    
    r = [0]*d
    for i in range(d):
        r[i] = r1[i]/r2[i] - 1
        med[i] += r[i] 
        medg[i] += r1[i]
        medr[i] += r2[i]
    print(c)

for i in range(d):
    med[i] = med[i]/count 
    medg[i] = medg[i]/count
    medr[i] = medr[i]/count
    
print(medg)
print(medr)
print(med)

# displaying the title 
label = 'Line Chart of the Correlations by Distance in the SIR Model'
plt.title(label, fontweight=10, pad='2.0') 
x = ('1', '2', '3', '4', '5', '6')
plt.plot(x, list([med[0], med[1], med[2], med[3], med[4], med[5]]))
# Create names for axis
plt.xlabel('Distances')
plt.ylabel('Correlation')
# Show plot
plt.figure()
plt.show()

# displaying the title 
label = 'Bar Plot of the Correlations by Distance in the SIR Model'
plt.title(label, fontweight=10, pad='2.0') 
y_plot = (0, med[0], med[1], med[2], med[3], med[4], med[5])
y_plot = pd.Series(y_plot)
plot1 = y_plot.plot.bar(grid=True, color='#607c8e')
plt.xlim([0.5, 6.5])
# Create names for axis
plt.xlabel('Distances')
plt.ylabel('Correlation')
# Show plot
plt.figure()
plt.show()

