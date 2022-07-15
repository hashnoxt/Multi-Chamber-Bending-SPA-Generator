"""
Implements a collection of gene classes

Genes support the following python operators:
    - + - calculates the phenotype resulting from the
      combination of a pair of genes

These genes work via classical Mendelian genetics
"""

import sys
from random import random, uniform

class BaseGene(object):
    """
    Base class from which all the gene classes are derived.

    You cannot use this class directly, because there are
    some methods that must be overridden.    
    """
    # each gene should have an object in
    # which its genotype should be stored
    value = None
    
    # probability of a mutation occurring
    mutProb = 0.01
    
    def __init__(self, value=None):
    
        # if value is not provided, it will be
        # randomly generated
        if value == None:
            value = self.randomValue()
    
        self.value = value
    
    def copy(self):
        """
        returns clone of this gene
        """
        return self.__class__(self.value)
    
    def __add__(self, other):
        """
        Combines two genes in a gene pair, to produce an effect
    
        This is used to determine the gene's phenotype
    
        This default method computes the arithmetic mean
        of the two genes.
    
        Override as needed
    
        Must be overridden
        """
        raise Exception("Method __add__ must be overridden")
    
    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.value)
    
    # def __cmp__(this, other):
    #     return cmp(this.value, other.value)
    # 
    def maybeMutate(self):
        if random() < self.mutProb:
            self.mutate()
    
    def mutate(self):
        """
        Perform a mutation on the gene
        
        You MUST override this in subclasses
        """
        raise Exception("method 'mutate' not implemented")
    
    def randomValue(self):
        """
        Generates a plausible random value
        for this gene.
        
        Must be overridden
        """
        raise Exception("Method 'randomValue' not implemented")
    
    
    def copy(self):
        """
        Produce an exact copy of this gene
    
        Shouldn't need to be overridden
        """
        return self.__class__(self.value)

class FloatGene(BaseGene):
    """
    A gene whose value is a floating point number

    Class variables to override:

        - mutAmt - default 0.1 - amount by which to mutate.
          The gene will will move this proportion towards
          its permissible extreme values

        - randMin - default -1.0 - minimum possible value
          for this gene. Mutation will never allow the gene's
          value to be less than this

        - randMax - default 1.0 - maximum possible value
          for this gene. Mutation will never allow the gene's
          value to be greater than this
    """
    # amount by which to mutate, will change value
    # by up to +/- this amount
    mutAmt = 0.1
    
    # used for random gene creation
    # override in subclasses
    randMin = -1.0
    randMax = 1.0
    
    def __add__(self, other):
        """
        Combines two genes in a gene pair, to produce an effect
    
        This is used to determine the gene's phenotype
    
        This class computes the arithmetic mean
        of the two genes' values, so is akin to incomplete
        dominance.
    
        Override if desired
        """
        return (self.value + other.value) / 2
    
    def mutate(self):
        """
        Mutate this gene's value by a random amount
        within the range, which is determined by
        multiplying self.mutAmt by the distance of the
        gene's current value from either endpoint of legal values
        
        perform mutation IN-PLACE, ie don't return mutated copy
        """
        if random() < 0.5:
            # mutate downwards
            self.value -= uniform(0, self.mutAmt * (self.value-self.randMin))
        else:
            # mutate upwards:
            self.value += uniform(0, self.mutAmt * (self.randMax-self.value))
    
    def randomValue(self):
        """
        Generates a plausible random value
        for this gene.
        
        Override as needed
        """
        min = self.randMin
        range = self.randMax - min
    
        return random() * range + min
    

class FloatGeneRandom(FloatGene):
    """
    Variant of FloatGene where mutation always randomises the value
    """
    def mutate(self):
        """
        Randomise the gene
    
        perform mutation IN-PLACE, ie don't return mutated copy
        """
        self.value = self.randomValue()
    

class FloatGeneMax(FloatGene):
    """
    phenotype of this gene is the greater of the values
    in the gene pair
    """
    def __add__(self, other):
        """
        produces phenotype of gene pair, as the greater of this
        and the other gene's values
        """
        return max(self.value, other.value)