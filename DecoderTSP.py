from brkga_mp_ipr.types import BaseChromosome
import numpy as np

class TSPDecoder:
    def __init__(self, numVertices, matrizCustos, prazos, arestas):
        self.numVertices = numVertices
        self.matrizCustos = matrizCustos
        self.prazos = prazos
        self.arestas = set(arestas)  
        self.M = 100000  # Valor grande para penalizar subrotas

    def decode(self, chromosome: BaseChromosome, rewrite: bool = False) -> float:
        # # Força a rota começar e terminar em 0
        permutation = [0] + sorted(range(1, self.numVertices), key=lambda k: chromosome[k-1]) + [0]

        total_cost = 0
        time = 0
        t = [0] * self.numVertices  # Cria a lista dos tempos de chegada de chegada

        visited = set()
        visited.add(0)  

        for k in range(len(permutation) - 1):
            u = permutation[k]
            v = permutation[k + 1]

            # aresta existe
            if (u, v) not in self.arestas:
                return float('inf')  # Penalidade absurda

            # tempo de chegada ao vértice v
            time += self.matrizCustos[u][v]
            t[v] = time

            #  prazo
            if v != 0 and t[v] > self.prazos[v]:
                return float('inf')  # Penalidade absurda
            
            # Restrições MTZ para eliminação de subrotas
            if t[u] + self.matrizCustos[u][v] - t[v] > self.M * (1 - int((u, v) in self.arestas)):
                return float('inf')  


            total_cost += self.matrizCustos[u][v]
            visited.add(v)

        # Garantia que a rota visite todos os vértices exatamente uma vez e retorne a 0
        if len(visited) != self.numVertices and permutation[0] != 0 and permutation[-1] != 0:
            return float('inf')
        

        return total_cost

