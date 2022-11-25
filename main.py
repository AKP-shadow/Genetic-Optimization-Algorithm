import math
import random
from collections import OrderedDict
import matplotlib.pyplot as plt


range_of_genes = (-600,600)
no_of_allele = 8
no_of_genes_per_chromosome = 5
no_of_chromosomes =6
no_of_parents_selected = 3


def RastigrinfitnessFunc( chromosome):
    fitness = 10*len(chromosome)
    mapping = MappingOfGenes()
 
    for i in range(len(chromosome)):
        fitness += mapping[chromosome[i]]**2 - (10*math.cos(2*math.pi*mapping[chromosome[i]]))
    return fitness

def GriewankfitnessFunc(chromosome):
    sum=0
    product=1
    mapping = MappingOfGenes()
    for i in range(len(chromosome)):
        sum+=mapping[chromosome[i]]**2
        product*=math.cos(mapping[chromosome[i]]/math.sqrt(i+1))
    return sum/4000 - product + 1

def Crossover(X,Y,no_of_allele):
    
    # RNG chromosome cutter
    cut = random.randint(1,no_of_allele*len(X)-1)
    if cut%no_of_allele==0:
        return X[:cut//no_of_allele]+Y[cut//no_of_allele:]
    else:
        gene_X = bin(X[cut//no_of_allele])[2:]
        gene_Y= bin(Y[cut//no_of_allele])[2:]
        while len(gene_X)!=no_of_allele:
            gene_X = '0'+gene_X
        while len(gene_Y)!=no_of_allele:
            gene_Y = '0'+gene_Y
            
        # Genes at CrossOver - Gene X and Gene Y
        # print(gene_X)
        # print(gene_Y)
        
        #Recombinant Gene
        # print(gene_X[:cut%no_of_allele])
        # print(gene_Y[cut%no_of_allele:])

        # Child Chromosome
        return X[:cut//no_of_allele] +[int(gene_X[:cut%no_of_allele]+gene_Y[cut%no_of_allele:],2)] + Y[cut//no_of_allele+1:]
    
    
def RankingOfSelectedParents(no_of_parents_selected,parents):
        
    # total_desired_fitness = 0
    # for parent in parents:
    #     total_desired_fitness+=fitness[parent]
    # print(total_desired_fitness)
    # probabilities={}
    
    rank={}
    total = no_of_parents_selected*(no_of_parents_selected+1)/2
    for i in range(no_of_parents_selected):
        rank[parents[i]] = (no_of_parents_selected - i)/total
        # probabilities[parents[i]] = fitness[parents[i]]/total_desired_fitness
    # print('probabilities:',probabilities)
    return rank
    
    
def SingleMutation(chromosome,no_of_allele,MutationRate=0.2):
    MutationAt = random.randint(0,len(chromosome)*no_of_allele-1)
    gene = bin(chromosome[(MutationAt)//no_of_allele])[2:]
    while len(gene)!=no_of_allele:
        gene = '0'+gene
    flip = str(random.choices([0,1],weights=[1-MutationRate,MutationRate],k=1)[0])
    gene = gene[:MutationAt%no_of_allele] + flip + gene[MutationAt%no_of_allele+1:]
    return chromosome[:MutationAt//no_of_allele] + [int(gene,2)] + chromosome[MutationAt//no_of_allele+1:]
        
def MappingOfGenes():
    mapping=[]
    s=range_of_genes[0]
    step = (range_of_genes[1]-range_of_genes[0])/(2**no_of_allele-1)
    for i in range(2**no_of_allele-1):
        mapping.append(int(s))
        s+=step
    mapping.append(int(s))   
    return mapping


def GetFitness(chromosomes,fitnessFunc):
    fitness = OrderedDict()
    if fitnessFunc == "Rastigrin":
        for i in range(len(chromosomes)):
            fitness[i] = RastigrinfitnessFunc(chromosomes[i])
        return fitness
    else:
        for i in range(len(chromosomes)):
            fitness[i] = GriewankfitnessFunc(chromosomes[i])
        return fitness      
    

def GA(chromosomes,n,fitnessFunc="Rastigrin",minimize=1):

    while n:
        fitness = GetFitness(chromosomes,fitnessFunc)
        if minimize:
            dominant_parents = sorted(fitness,key=lambda x:fitness[x])
        else:
            dominant_parents = sorted(fitness,key=lambda x:fitness[x],reverse=True)
        recessive_parents = dominant_parents[no_of_parents_selected:]
        dominant_parents = dominant_parents[:no_of_parents_selected]
 
        rank = RankingOfSelectedParents(no_of_parents_selected,dominant_parents)
        naturally_selected_parent_pairs=[]
        k=0
        while k<no_of_parents_selected:
            parent_s = random.choices(dominant_parents,weights=rank.values(),k=2)
            
            if parent_s[0]!=parent_s[1]:
                k+=1
                naturally_selected_parent_pairs.append(parent_s)
            
        # print(naturally_selected_parent_pairs)
        
        # Crossover of parents
        children=[]
        for i in range(no_of_parents_selected):
            child = Crossover(chromosomes[naturally_selected_parent_pairs[i][0]],chromosomes[naturally_selected_parent_pairs[i][1]],no_of_allele)
            # Single Mutation
            children.append(SingleMutation(child,no_of_allele))

        for ele in sorted(recessive_parents,reverse=True):
            del chromosomes[ele]
        
        chromosomes+=children
        n-=1
        
    return min(GetFitness(chromosomes,fitnessFunc).values())




initial_chromosomes=[]
#random Chromosome selection
for i in range(no_of_chromosomes):
    c=[]
    for j in range(no_of_genes_per_chromosome):
        c.append(random.randint(0,2**no_of_allele-1))
    initial_chromosomes.append(c)
    
print("Initial group of chromosomes")
for i in range(len(initial_chromosomes)):
    print(f"Chromosome#{i+1}: ",initial_chromosomes[i])
print("Range of Genes: ",range_of_genes)
print("No of allele: ",no_of_allele)
print("No of chromosomes: ",no_of_chromosomes)
print("No of Genes per chromosome: ",no_of_genes_per_chromosome)
print("Percentage of parents selected: ",no_of_parents_selected/no_of_chromosomes*100,"%")
generations = 200
print("No of Generations: ",generations)
for fitnessFunc in ['Rastigrin','Griewank']:
    for j in range(2):
        recombinant_min_values=[]
        for i in range(generations):
            recombinant_min_values.append(GA(initial_chromosomes,i,fitnessFunc,minimize=j))
    
        print(f"{fitnessFunc} - {('Maximum','Minimum')[j]} value: ",recombinant_min_values[-1])
        plt.plot(range(generations),recombinant_min_values)
        plt.title(f"{fitnessFunc} - {('Maximize','Minimize')[j]}")
        plt.xlabel("'x' generations")
        plt.ylabel("Selected optmimum value after 'x' generations")
        # plt.show()
        plt.savefig(f'GAplt_{fitnessFunc}_{j}.jpeg')
        plt.close()

