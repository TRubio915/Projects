# valueIterationAgents.py
# -----------------------
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


# valueIterationAgents.py
# -----------------------
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


import mdp, util

from learningAgents import ValueEstimationAgent
import collections

class ValueIterationAgent(ValueEstimationAgent):
    """
        * Please read learningAgents.py before reading this.*
        A ValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs value iteration
        for a given number of iterations using the supplied
        discount factor.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 100):
        """
          Your value iteration agent should take an mdp on
          construction, run the indicated number of iterations
          and then act according to the resulting policy.
          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state, action, nextState)
              mdp.isTerminal(state)
        """
        self.mdp = mdp
        self.discount = discount
        self.iterations = iterations
        self.values = util.Counter() # A Counter is a dict with default 0
        self.runValueIteration()
    def runValueIteration(self):
        # Write value iteration code here
        "*** YOUR CODE HERE ***"
        value = util.Counter()
        states = self.mdp.getStates()
        #loop through each iterations
        #copy values over
        for it in range(0,self.iterations):
            value = self.values.copy()
            #for each state get the actions
            for state in states:
                actions = self.mdp.getPossibleActions(state)
                myvals = []
                #if the state is terminal then set it to 0
                if self.mdp.isTerminal(state):
                    self.values[state] = 0
                #otherwise go through actions
                else:
                    #for each action get the transition states and probabilities
                    for action in actions:
                        TSAPs = self.mdp.getTransitionStatesAndProbs(state, action)
                        newval = 0
                        #for each transition state and probability, calculate the value
                        for TSAP in TSAPs:
                            newval = (newval + TSAP[1] * (self.mdp.getReward(state, action, TSAP[0]) + self.discount * value[TSAP[0]]))
                        myvals.append(newval)
                    self.values[state] = max(myvals)

    def getValue(self, state):
        """
          Return the value of the state (computed in __init__).
        """
        return self.values[state]

    def computeQValueFromValues(self, state, action):
        """
          Compute the Q-value of action in state from the
          value function stored in self.values.
        """
        "*** YOUR CODE HERE ***"
        value = 0
        TSAPs = self.mdp.getTransitionStatesAndProbs(state, action)
        #loop through transition states and probabiliies
        for TSAP in TSAPs:
            #calculate a new total qvalue from the values
             value = value + (TSAP[1] * (self.mdp.getReward(state, action, TSAP[0]) + self.discount * self.values[TSAP[0]]))
        return value
        util.raiseNotDefined()

    def computeActionFromValues(self, state):
        """
          The policy is the best action in the given state
          according to the values currently stored in self.values.
          You may break ties any way you see fit.  Note that if
          there are no legal actions, which is the case at the
          terminal state, you should return None.
        """
        "*** YOUR CODE HERE ***"
        #if the state is the termainl, return none because their is no action
        if self.mdp.isTerminal(state):
            return None
        else:
            max_action = 0
            max_value = float("-inf")
            actions = self.mdp.getPossibleActions(state)
            #for each action want to loop through transition states and probabilities
            for action in actions:
                value = 0
                TSAPs = self.mdp.getTransitionStatesAndProbs(state, action)
                #loop through each transition state and probability
                #calculate new value per transition state and probability
                for TSAP in TSAPs:
                    value = value + (TSAP[1] * (self.mdp.getReward(state, action, TSAP[0]) + self.discount * self.values[TSAP[0]]))
                #if current value is greater than the max, update it
                if value > max_value:
                    max_action = action
                    max_value = value
            return max_action

    def getPolicy(self, state):
        return self.computeActionFromValues(state)

    def getAction(self, state):
        "Returns the policy at the state (no exploration)."
        return self.computeActionFromValues(state)

    def getQValue(self, state, action):
        return self.computeQValueFromValues(state, action)

class AsynchronousValueIterationAgent(ValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        An AsynchronousValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs cyclic value iteration
        for a given number of iterations using the supplied
        discount factor.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 1000):
        """
          Your cyclic value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy. Each iteration
          updates the value of only one state, which cycles through
          the states list. If the chosen state is terminal, nothing
          happens in that iteration.

          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state)
              mdp.isTerminal(state)
        """
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        #Looping through iterations
        for i in range(self.iterations):
            #Variable assignment
            curState = self.mdp.getStates()[i %  len(self.mdp.getStates())] #Current states from MDP
            topVal = self.computeActionFromValues(curState) #Best value found
            #check if our best value is none, if so assign our current value '0'
            if topVal is None:
                curVal = 0
            #If not compute and assign our current value
            else:
                curVal = self.computeQValueFromValues(curState, topVal)
            #Assign value of our current state to computed value
            self.values[curState] = curVal

    #Function to return a counter variable containing our qValues
    def qValueFind(self, state):
        #Declarations
        actList = self.mdp.getPossibleActions(state) #A list of our state's possible actions
        values = util.Counter() #Counter to hold qValues
        #Begin loop to fill counter with values computed from state and action
        for act in actList:
            values[act] = self.computeQValueFromValues(state, act)
        #return counter
        return values


class PrioritizedSweepingValueIterationAgent(AsynchronousValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        A PrioritizedSweepingValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs prioritized sweeping value iteration
        for a given number of iterations using the supplied parameters.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 100, theta = 1e-5):
        """
          Your prioritized sweeping value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy.
        """
        self.theta = theta
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        "*** YOUR CODE HERE ***"

