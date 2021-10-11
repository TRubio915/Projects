# qlearningAgents.py
# ------------------
# Licensing Information:  You are free to use or extend these projects for
# educational purposes provided that (1) you do not distribute or publish
# solutions, (2) you retain this notice, and (3) you provide clear
# attribution to UC Berkeley, including a link to http://ai.berkeley.edu.
# 
# Attribution Information: The Pacman AI projects were developed at UC Berkeley.
# The core projects and autograders were primarily created by John DeNero
# (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and
# Pieter Abbeel (pabbeel@cs.berkeley.edu).


from game import *
from learningAgents import ReinforcementAgent
from featureExtractors import *

import random,util,math

class QLearningAgent(ReinforcementAgent):
    """
      Q-Learning Agent

      Functions you should fill in:
        - computeValueFromQValues
        - computeActionFromQValues
        - getQValue
        - getAction
        - update

      Instance variables you have access to
        - self.epsilon (exploration prob)
        - self.alpha (learning rate)
        - self.discount (discount rate)

      Functions you should use
        - self.getLegalActions(state)
          which returns legal actions for a state
    """
    def __init__(self, **args):
        "You can initialize Q-values here..."
        ReinforcementAgent.__init__(self, **args)
        self.qvalues = util.Counter() #Initializing qValues using a Counter utility

    def getQValue(self, state, action):
        """
          Returns Q(state,action)
          Should return 0.0 if we have never seen a state
          or the Q node value otherwise
        """
        #Return calculated qValue for provided state and action
        return self.qvalues[(state, action)]


    def computeValueFromQValues(self, state):
        """
          Returns max_action Q(state,action)
          where the max is over legal actions.  Note that if
          there are no legal actions, which is the case at the
          terminal state, you should return a value of 0.0.
        """
        #Declarations
        val = float('-inf') #initialize value to negative infinity
        actList = self.getLegalActions(state) #List of legal actions from state
        # Begin looping through each legal action
        for act in actList:
            #If qValue is greater than or equal to current value, replace
            if self.getQValue(state, act) >= val:
                val = self.getQValue(state, act)
        #Return value as long as it is not equal to negative infinity (unchanged)
        return val if val != (float('-inf')) else 0.0

    def computeActionFromQValues(self, state):
        """
          Compute the best action to take in a state.  Note that if there
          are no legal actions, which is the case at the terminal state,
          you should return None.
        """
        #Declarations
        val = float('-inf') #initialize value to negative infinity
        actList = self.getLegalActions(state) #List of legal actions from state
        #Begin looping through each legal action
        for act in actList:
            #if qValue is greater than current value
            if self.getQValue(state, act) > val:
                #Assign qValue to current value
                val = self.getQValue(state, act)
                topAct = act #We have found a new "best" action
        list = [] #Create empty list
        #Second loop through each legal action
        for act in actList:
            #if qValue is equal to our current value, add to list
            if self.getQValue(state, act) == val:
                list.append(act)
        #if our current value does not equal negative inf, AKA has been changed
        if val != float('-inf'):
            #Assign best action a random action from list
            topAct = random.choice(list)
        #Else change to NULL
        else:
            topAct = None
        #return best action found
        return topAct

    def getAction(self, state):
        """
          Compute the action to take in the current state.  With
          probability self.epsilon, we should take a random action and
          take the best policy action otherwise.  Note that if there are
          no legal actions, which is the case at the terminal state, you
          should choose None as the action.
          HINT: You might want to use util.flipCoin(prob)
          HINT: To pick randomly from a list, use random.choice(list)
        """
        # Pick Action
        legalActions = self.getLegalActions(state)
        action = None
        if util.flipCoin(self.epsilon):
            action = random.choice(legalActions)
        else:
            action = self.getPolicy(state)
        #util.raiseNotDefined()
        return action

    def update(self, state, action, nextState, reward):
        """
          The parent class calls this to observe a
          state = action => nextState and reward transition.
          You should do your Q-Value update here

          NOTE: You should never call this function,
          it will be called on your behalf
        """
        #Begin update of qVales
        alphaQ = (self.getQValue(state, action)) * (1 - self.alpha) #Get qValue from provided state and action, multiplied times 1 - alpha
        rewardQ = (reward + self.discount * self.computeValueFromQValues(nextState)) * (self.alpha) #Computing qValues with reward and discount
        self.qvalues[(state, action)] = alphaQ + rewardQ #Assign sum of two variables to update qValues

    def getPolicy(self, state):
        return self.computeActionFromQValues(state)

    def getValue(self, state):
        return self.computeValueFromQValues(state)


class PacmanQAgent(QLearningAgent):
    "Exactly the same as QLearningAgent, but with different default parameters"

    def __init__(self, epsilon=0.05,gamma=0.8,alpha=0.2, numTraining=0, **args):
        """
        These default parameters can be changed from the pacman.py command line.
        For example, to change the exploration rate, try:
            python pacman.py -p PacmanQLearningAgent -a epsilon=0.1

        alpha    - learning rate
        epsilon  - exploration rate
        gamma    - discount factor
        numTraining - number of training episodes, i.e. no learning after these many episodes
        """
        args['epsilon'] = epsilon
        args['gamma'] = gamma
        args['alpha'] = alpha
        args['numTraining'] = numTraining
        self.index = 0  # This is always Pacman
        QLearningAgent.__init__(self, **args)

    def getAction(self, state):
        """
        Simply calls the getAction method of QLearningAgent and then
        informs parent of action for Pacman.  Do not change or remove this
        method.
        """
        action = QLearningAgent.getAction(self,state)
        self.doAction(state,action)
        return action


class ApproximateQAgent(PacmanQAgent):
    """
       ApproximateQLearningAgent

       You should only have to overwrite getQValue
       and update.  All other QLearningAgent functions
       should work as is.
    """
    def __init__(self, extractor='IdentityExtractor', **args):
        self.featExtractor = util.lookup(extractor, globals())()
        PacmanQAgent.__init__(self, **args)
        self.weights = util.Counter()

    def getWeights(self):
        return self.weights

    def getQValue(self, state, action):
        """
          Should return Q(state,action) = w * featureVector
          where * is the dotProduct operator
        """
        #Declarations
        featList = self.featExtractor.getFeatures(state, action) #Vector of features for provided state and action
        weight = self.getWeights() #Assign weights
        val = weight * featList #qValue to be returned
        #Return calculated value
        return val

    def update(self, state, action, nextState, reward):
        """
           Should update your weights based on transition
        """
        #Declarations
        maxQ = self.computeValueFromQValues(nextState) #Max qValue taken from provide next state
        featList = self.featExtractor.getFeatures(state, action) #Vector of features for provided state and action
        actQVal = self.getQValue(state, action) #Q value for provided action
        #Begin looping through each feature in vector
        for feat in featList:
            #Update weights for feature
            update = ((self.discount * maxQ + reward) - actQVal) * self.alpha * featList[feat] #Calculated update
            self.weights[feat] += update #Add update

    # Update override
    def update(self, state, action, nextState, reward):
        """
           Should update your weights based on transition
        """
        #Declarations
        featList = self.featExtractor.getFeatures(state, action) #Vector of features for provided state and action
        actQVal = self.getQValue(state, action) #Q value for provided action
        difference = (reward + self.discount * self.computeValueFromQValues(nextState)) - actQVal #Holds calculated difference term
        #Begin looping through each feature in vector
        for feat in featList:
            #Update our weight incorporating difference
            update = featList[feat] * difference * self.alpha #Calculated update
            self.weights[feat] += update #Add update

    def final(self, state):
        "Called at the end of each game."
        # call the super-class final method
        PacmanQAgent.final(self, state)

        # did we finish training?
        if self.episodesSoFar == self.numTraining:
            # you might want to print your weights here for debugging
            "*** YOUR CODE HERE ***"
            pass


