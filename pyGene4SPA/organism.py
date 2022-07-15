"""
Implements classes for entire organisms

Organisms can be mated by the '+' operator, which
produces a child organism.

Subclasses of Organism must override the following methods:
    - fitness - returns a float value representing the
      organism's fitness - a value from 0.0 to infinity, where
      lower is better
"""
import math
from random import random, choice
from gene import BaseGene

class Organism(object):
    """
    Simple genetic algorithms organism

    Contains only single genes, not pairs (ie, single-helix)
    
    Note - all organisms are hermaphrodites, which
    can reproduce by mating with another.
    In this implementation, there is no gender.

    Class variables (to override) are:
        - genome - a dict mapping gene names to gene classes

        - mutateOneOnly - default False - dictates whether mutation
          affects one randomly chosen gene unconditionally, or
          all genes subject to the genes' individual mutation settings

        - crossoverRate - default .5 - proportion of genes to
          split out to first child in each pair resulting from
          a mating

    Python operators supported:
        - + - mates two organism instances together, producing a child
        - [] - returns the value of the gene of a given name
        - <, <=, >, >= - compares fitness value to that of another instance
    """
    # dict which maps genotype names to gene classes
    genome = {}
    
    #default value type for a gene
    geneVType = type(2.0)
    # dictates whether mutation affects one randomly chosen
    # gene unconditionally, or all genes subject to the genes'
    # own mutation settings
    
    mutateOneOnly = False
    
    # proportion of genes to split out to first child
    crossoverRate = 0.5
    
    def __init__(self, *arg, **kw):
        """
        Initialises this organism randomly,
        or from a set of named gene keywords
    
        Arguments:
            - gamete1, gamete2 - a pair of gametes from which
              to take the genes comprising the new organism.
              May be omitted.
        
        Keywords:
            - keyword names are gene names within the organism's
              genome, and values are either:
                  - instances of a Gene subclass, or
                  - a Gene subclass (in which case the class will
                    be instantiated to form a random gene object)
    
        Any gene names in the genome, which aren't given in the
        constructor keywords, will be added as random instances
        of the respective gene class. (Recall that all Gene subclasses
        can be instantiated with no arguments to create a random
        valued gene).
        """
        # the set of genes which comprise this organism
        self.genes = {}
    
        # remember the gene count
        self.numgenes = len(self.genome)
        if len(arg)<=self.numgenes:
            arg = list(arg) + [None]*(self.numgenes-len(arg))
        else:
            arg = arg[0:self.numgenes]
        # we're being fed a set of zero or more genes
        p = 0
        for name, cls in self.genome.items():
    
            # set genepair from given arg, or default to a
            # new random instance of the gene
            gene = kw.get(name, cls)
            flag = (type(arg[p])==self.geneVType)
            # if we're handed a gene class instead of a gene object
            # we need to instantiate the gene class
            # to form the needed gene object
            try:
                if issubclass(gene, BaseGene):
                    gene = (gene(arg[p]) if flag else gene())
                else:
                    gene = (cls(arg[p]) if flag else cls())
            except:
                pass
            # either way, we should have a valid gene now
            # add in the gene to our genotype
            self.genes[name] = gene
            p = p +1

        #record the fittness result for the org
        value = kw.get('FIT')
        if value:
            self.fitnessValue = value
        else:
            self.fitnessValue = self.fitness()
    
    def copy(self):
        """
        returns a deep copy of this organism
        """
        genes = {}
        for name, gene in self.genes.items():
            genes[name] = gene.copy()
        genes['FIT'] = self.fitnessValue
        return self.__class__(**genes)
    
    def mate(self, partner):
        """
        Mates this organism with another organism to
        produce two entirely new organisms via random choice
        of genes from this or the partner
        """
        genotype1 = {}
        genotype2 = {}
    
        # gene by gene, we assign our and partner's genes randomly
        for name, cls in self.genome.items():
            
            ourGene = self.genes.get(name, None)
            if not ourGene:
                ourGene = cls()
    
            partnerGene = self.genes.get(name, None)
            if not partnerGene:
                partnerGene = cls()
                
            # randomly assign genes to first or second child
            if random() < self.crossoverRate:
                genotype1[name] = ourGene
                genotype2[name] = partnerGene
            else:
                genotype1[name] = partnerGene
                genotype2[name] = ourGene
        
        # got the genotypes, now create the child organisms
        child1 = self.__class__(**genotype1)
        child2 = self.__class__(**genotype2)
        # done
        return (child1, child2)
    
    def __getitem__(self, item):
        """
        allows shorthand for querying the phenotype
        of this organism
        """
        return self.genes[item].value
    
    def phenotype(self, geneName=None):
        """
        Returns the phenotype resulting from a
        given gene, OR the total phenotype resulting
        from all the genes
        
        tries to invoke a child class' method
        called 'phen_<name>'
        """
        # if no gene name specified, build up an entire
        # phenotype dict
        if geneName == None:
            phenotype = {}
            for name, cls in self.genome.items():
                val = self.phenotype(name)
                if not phenotype.has_key(name):
                    phenotype[name] = []
                phenotype[name].append(val)
    
            # got the whole phenotype now
            return phenotype
    
        # just getting the phenotype for one gene pair
        return self.genes[geneName]
    
    def mutate(self):
        """
        Implement the mutation phase, invoking
        the stochastic mutation method on each
        component gene
        
        Does not affect this organism, but returns a mutated
        copy of it
        """
        mutant = self.copy()
        
        if self.mutateOneOnly:
            # unconditionally mutate just one gene
            gene = choice(mutant.genes.values())
            gene.mutate()
    
        else:
            # conditionally mutate all genes
            for gene in mutant.genes.values():
                gene.maybeMutate()
    
        return mutant
    
    def __add__(self, partner):
        """
        Allows '+' operator for sexual reproduction
    
        Returns a whole new organism object, whose
        gene pair for each gene name are taken as one
        gene randomly selected from each parent
        """
        return self.mate(partner)
    
    def fitness(self):
        """
        Return the fitness level of this organism, as a float
        
        Should return a number from 0.0 to infinity, where
        0.0 means 'perfect'
    
        Organisms should evolve such that 'fitness' converges
        to zero.
        
        This method must be overridden
        """
        raise Exception("Method 'fitness' not implemented")
    
    def duel(self, opponent):
        """
        Duels this organism against an opponent
        
        Returns -1 if this organism loses, 0 if it's
        a tie, or 1 if this organism wins
        """
        return cmp(self.fitnessValue, opponent.fitnessValue)
    
    def __cmp__(self, other):
        """
        Convenience method which invokes duel
        
        Allows lists of organisms to be sorted
        """
        return self.duel(other)
    
    def __repr__(self):
        """
        Delivers a minimal string representation
        of this organism.
        
        Override if needed
        """
        return "<%s:%s>" % (self.__class__.__name__, self.fitnessValue)