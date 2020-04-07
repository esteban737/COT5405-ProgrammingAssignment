import sys
import numpy as np
import random
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

t = 50000
q = lambda p : (1 - p)
tintervals = np.arange(1, t + 1, 1)
kintervals = 0
tick_spacing_time = t/5

class Network:
    def __init__ (self):
        self.current = 1
        self.n = 1.0
        self.m = 1
        self.nt = [0]
        self.mt = [0]
        self.nodesK = [1]
        self.graph = {1: [1]}

    def addNode(self, target):
        if self.n != 0:
            self.graph[self.current+1] = [target]
            targetL = len(self.graph[target])

            if targetL != 0:
                self.nodesK[targetL - 1] -= 1

            if targetL == len(self.nodesK):
                self.nodesK.append(1)
            else:
                self.nodesK[targetL] += 1

            self.graph[target].append(self.current+1)
        else:
            self.graph[target] = [target]

        self.n += 1
        self.m += 1
        self.nodesK[0] += 1
        self.current += 1

    def deleteNode(self, node):
        if node == 0:
            return

        targets = self.graph[node]

        if len(targets) != 0:
            self.nodesK[len(targets) - 1] -= 1
            
        for target in targets:
            if node != target:
                l = len(self.graph[target])
                if l - 1 == 0:
                    self.nodesK[0] -= 1
                    self.graph[target].remove(node)
                    self.m -= 1
                else:
                    self.nodesK[l - 1] -= 1
                    self.nodesK[l - 2] += 1
                    self.graph[target].remove(node)
                    self.m -= 1
            else:
                self.m -= 1
        del self.graph[node]
        self.n -= 1

    def updateStats(self):
        self.nt.append(self.n)
        self.mt.append(self.m)

        
class Numerical: 
    def __init__ (self):
        self.nt = [] * t
        self.mt = [] * t
        self.nodesk = []

    def numNodes(self, p, t):
        res = (p - q(p)) * t + 2 * q(p)
        self.nt.append(res)

    def numEdges(self, p, t):
        res = p * t * (p - q(p))
        self.mt.append(res)
    
    def dist(self, p, k, n):
        res = pow(k, (-1) - (2 * p / (2 * p - 1)))
        self.nodesk.append(res)

def main():
    setSimNodes = []
    setSimEdges = []

    setNumNodes = []
    setNumEdges = []

    numDistModel = []
    numDist = []

    pVals = [0.6, 0.75, 0.9]
    kintervals = []
    reps = 1

    for p in pVals:
        results = testingCycle(p, reps)
        setSimNodes.append(results[0][0])
        setSimEdges.append(results[0][1])
        setNumNodes.append(results[1][0])
        setNumEdges.append(results [1][1])

    numDistModel, kintervals, net = runDistSimulation()
    numDist = runNumerical(kintervals, net.n)

    graphVars(setSimNodes, setNumNodes, 'Nodes')
    graphVars(setSimEdges, setNumEdges, 'Edges')
    graphDist(numDistModel, numDist, kintervals)

def updateNetwork(p, net):
    choice = random.random()
    if choice <= (p): 
        birth(net)
    else: 
        death(net)
   
def testingCycle(p, reps):
    numResults = []
    simResults = []

    networks = map(lambda x: Network(), range(1, reps + 1))
    netNodes = [0]
    netEdges = [0]
    num = Numerical()

    for t in tintervals:
        num.numNodes(p, t)
        num.numEdges(p, t)
        map(lambda net: multiSim(t, p, net), networks)
    numResults.append(num.nt)
    numResults.append(num.mt)
    for i in range(1, 6):
        sumNodes = 0
        sumEdges = 0
        for net in networks:
            sumNodes += net.nt[i]
            sumEdges += net.mt[i]
        netNodes.append(sumNodes/reps)
        netEdges.append(sumEdges/reps)


    simResults.append(netNodes)
    simResults.append(netEdges)
    return [simResults, numResults]


def multiSim(t, p, net):
    updateNetwork(p,net)
    if t % tick_spacing_time == 0:
                net.updateStats()


def runDistSimulation():
    kintervals = []
    normalizedData = []
    net = Network()
    for t in tintervals:
        updateNetwork(0.8, net)
    for k in net.nodesK:
        res = k/net.n
        if res != 0:
            normalizedData.append(k/net.n)
    kintervals = np.arange(1, len(normalizedData) + 1, 1)
    return [normalizedData, kintervals, net]

def runNumerical(kintervals, n):
    num = Numerical()
    for k in kintervals:
        num.dist(0.8, k, n)
    return num.nodesk

def birth(net):
    nodes = []
    probs  = []

    for key, val in net.graph.items():
        nodes.append(key)
        
        if net.m == 0:
            probs.append(1/len(net.graph))
        else:
            probs.append(len(val) / (2 * net.m))

    if len(nodes) != 0:
        target = np.random.choice(nodes, 1, probs)[0]
        net.addNode(target)
    else:
        net.addNode(net.current + 1)


def death(net):
    nodes = []
    probs  = []

    for key, val in net.graph.items():
        nodes.append(key)
        if net.n == 2 and net.m == 2:
            probs.append(0.5)
        elif net.m == 0:
            probs.append(1/len(net.graph))
        else:
            probs.append((net.n - len(val)) / (pow(net.n, 2) - (2 * net.m)))

    if len(nodes) != 0:
        node = np.random.choice(nodes, 1, probs)[0]
        net.deleteNode(node)
    else:
        net.deleteNode(0)
    
def graphVars(network, num, data):
    netData = np.linspace(0, t, num = 6)
    fig, ax = plt.subplots(1,1)
    
    ax.plot(netData, network[0], 'ro', label='p = 0.6')
    ax.plot(tintervals, num[0], 'r')
    ax.plot(netData, network[1], 'bs', label='p = 0.75')
    ax.plot(tintervals, num[1], 'b')
    ax.plot(netData, network[2], 'g^', label='p = 0.9')
    ax.plot(tintervals, num[2], 'g')
    plt.ylabel(data)
    plt.xlabel("t")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_time))
    plt.show()

def graphDist(network, num, kintervals):
    fig, ax = plt.subplots(1,1)
    ax.loglog(kintervals, network, 'r--', label='Model')
    ax.loglog(kintervals, num, 'b', label='Analytical')
    plt.ylabel("P(k)")
    plt.xlabel("k")
    ax.legend()
    plt.show()

main()
 
