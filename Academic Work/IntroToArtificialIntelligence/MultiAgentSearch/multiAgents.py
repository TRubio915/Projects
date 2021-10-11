# multiAgents.py
# --------------
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


from util import manhattanDistance
from game import Directions
import random, util

from game import Agent

class ReflexAgent(Agent):
    """
    A reflex agent chooses an action at each choice point by examining
    its alternatives via a state evaluation function.

    The code below is provided as a guide.  You are welcome to change
    it in any way you see fit, so long as you don't touch our method
    headers.
    """


    def getAction(self, gameState):
        """
        You do not need to change this method, but you're welcome to.

        getAction chooses among the best options according to the evaluation function.

        Just like in the previous project, getAction takes a GameState and returns
        some Directions.X for some X in the set {NORTH, SOUTH, WEST, EAST, STOP}
        """
        # Collect legal moves and successor states
        legalMoves = gameState.getLegalActions()

        # Choose one of the best actions
        scores = [self.evaluationFunction(gameState, action) for action in legalMoves]
        bestScore = max(scores)
        bestIndices = [index for index in range(len(scores)) if scores[index] == bestScore]
        chosenIndex = random.choice(bestIndices) # Pick randomly among the best

        "Add more of your code here if you want to"

        return legalMoves[chosenIndex]

    def evaluationFunction(self, currentGameState, action):
        """
        Design a better evaluation function here.

        The evaluation function takes in the current and proposed successor
        GameStates (pacman.py) and returns a number, where higher numbers are better.

        The code below extracts some useful information from the state, like the
        remaining food (newFood) and Pacman position after moving (newPos).
        newScaredTimes holds the number of moves that each ghost will remain
        scared because of Pacman having eaten a power pellet.

        Print out these variables to see what you're getting, then combine them
        to create a masterful evaluation function.
        """
        # Useful information you can extract from a GameState (pacman.py)
        successorGameState = currentGameState.generatePacmanSuccessor(action)
        newPos = successorGameState.getPacmanPosition()
        newFood = successorGameState.getFood()
        newGhostStates = successorGameState.getGhostStates()
        newScaredTimes = [ghostState.scaredTimer for ghostState in newGhostStates]

        "*** YOUR CODE HERE ***"
        #ghost code
        ghostdis = 2
        howcloseistheghost = 0
        #loop through each ghost to update where it is located
        for ghost in newGhostStates:
            distance = util.manhattanDistance(newPos, ghost.getPosition())
            #if the ghost is either distance 1 or 0, update howcloseistheghost for the fact that it is closer to pacman
            if distance <= 1:
                howcloseistheghost = howcloseistheghost + 1
            #update the specific ghost disatnae with the distance
            ghostdis = ghostdis + distance
        #determine part of the score with ghostdistance and how close the ghost is
        ghostscore = ((1/ghostdis)+howcloseistheghost)
        #food code
        #get newfood into a list of tuples
        lstNewFood = newFood.asList()
        #give food distance a starting of -1 since the place is filled with food
        fooddis = -1
        #loop through food pieces in the list of food locations
        for foodpiece in lstNewFood:
            distance = util.manhattanDistance(newPos, foodpiece)
            #if the food is one away or greater than/equal to the distance of pacmans position and food, then update it
            if fooddis == -1 or fooddis >= distance:
                fooddis = distance
        #determine part of the score with one over the distance of the food
        foodscore = (1/fooddis)
        #determine final score
        final_score = successorGameState.getScore() + foodscore - ghostscore
        return final_score

def scoreEvaluationFunction(currentGameState):
    """
    This default evaluation function just returns the score of the state.
    The score is the same one displayed in the Pacman GUI.

    This evaluation function is meant for use with adversarial search agents
    (not reflex agents).
    """
    return currentGameState.getScore()

class MultiAgentSearchAgent(Agent):
    """
    This class provides some common elements to all of your
    multi-agent searchers.  Any methods defined here will be available
    to the MinimaxPacmanAgent, AlphaBetaPacmanAgent & ExpectimaxPacmanAgent.

    You *do not* need to make any changes here, but you can if you want to
    add functionality to all your adversarial search agents.  Please do not
    remove anything, however.

    Note: this is an abstract class: one that should not be instantiated.  It's
    only partially specified, and designed to be extended.  Agent (game.py)
    is another abstract class.
    """

    def __init__(self, evalFn = 'scoreEvaluationFunction', depth = '2'):
        self.index = 0 # Pacman is always agent index 0
        self.evaluationFunction = util.lookup(evalFn, globals())
        self.depth = int(depth)

class MinimaxAgent(MultiAgentSearchAgent):
    """
    Your minimax agent (question 2)
    """

    def getAction(self, gameState):
        """
        Returns the minimax action from the current gameState using self.depth
        and self.evaluationFunction.

        Here are some method calls that might be useful when implementing minimax.

        gameState.getLegalActions(agentIndex):
        Returns a list of legal actions for an agent
        agentIndex=0 means Pacman, ghosts are >= 1

        gameState.generateSuccessor(agentIndex, action):
        Returns the successor game state after an agent takes an action

        gameState.getNumAgents():
        Returns the total number of agents in the game

        gameState.isWin():
        Returns whether or not the game state is a winning state

        gameState.isLose():
        Returns whether or not the game state is a losing state
        """
        def minimaxDispatch(gameState, depth):
            #Declarations
            actions = gameState.getLegalActions() #variable for state's legal actions
            topAction = Directions.STOP #Initialize best action to neutral STOP
            topScore = float("-inf") #Initialize best score to -infinity
            curScore = topScore #Holds our current score, initialized equal to top score

            #Base Case
            if gameState.isWin() or gameState.isLose():
                return gameState.getScore()
            #Begin loop
            for action in actions:
                curScore = maxValue(gameState.generateSuccessor(0, action), depth, 1)
                #if our current score surpasses our best score, replace best score and best action
                if curScore > topScore:
                    topAction = action
                    topScore = curScore
            #If our depth is 0, return top action found
            if depth == 0:
                return topAction
            #Otherwise return best score found
            else:
                return topScore

        def maxValue(gameState, depth, ghost):
            #Declarations
            nxtGhost = ghost + 1 #Variable to iterate through ghosts during their turn
            actions = gameState.getLegalActions(ghost) #variable for state's legal actions
            topScore = float("inf") #Initialize best score to +infinity
            curScore = topScore #Holds our current score, initialized equal to top score
            #Base Cases
            if gameState.isWin() or gameState.isLose():
                return gameState.getScore()
            if ghost == gameState.getNumAgents()-1:
                nxtGhost = 0
            #Begin Loop
            for action in actions:
                #If last ghost in ghosts' turn
                if nxtGhost == 0:
                    if depth == self.depth - 1:
                        curScore = self.evaluationFunction(gameState.generateSuccessor(ghost, action))
                    else:
                        curScore = minimaxDispatch(gameState.generateSuccessor(ghost, action), depth+1)
                #otherwise we recursively call our max function
                else:
                    curScore = maxValue(gameState.generateSuccessor(ghost, action), depth, nxtGhost)
                #if our current score is higher than our best score, replace
                if curScore < topScore:
                    topScore = curScore
            #Return current best score
            return topScore
        #Recursively call our dispatch function
        return minimaxDispatch(gameState,0)

class AlphaBetaAgent(MultiAgentSearchAgent):
    """
    Your minimax agent with alpha-beta pruning (question 3)
    """

    def getAction(self, gameState):
        """
        Returns the minimax action using self.depth and self.evaluationFunction
        """
        "*** YOUR CODE HERE ***"
        def minimaxDispatch(gameState, depth, alpha, beta):
            #Declarations
            actions = gameState.getLegalActions() #variable for state's legal actions
            topAction = Directions.STOP #Initialize best action to neutral STOP
            topScore = float("-inf") #Initialize best score to -infinity
            curScore = topScore #Holds our current score, initialized equal to top score
            #Base Case
            if gameState.isWin() or gameState.isLose():
                return gameState.getScore()
            #Begin loop
            for action in actions:
                curScore = maxValue(gameState.generateSuccessor(0, action), depth, 1, alpha, beta)
                #if our current score surpasses our best score, replace best score and best action
                if curScore > topScore:
                    topAction = action
                    topScore = curScore
                #update alpha while topscore changes
                alpha = max(alpha, topScore)
                if topScore > beta:
                    return topScore
            #If our depth is 0, return top action found
            if depth == 0:
                return topAction
            #Otherwise return best score found
            else:
                return topScore
        def maxValue(gameState, depth, ghost, alpha, beta):
            #Declarations
            nxtGhost = ghost + 1 #Variable to iterate through ghosts during their turn
            actions = gameState.getLegalActions(ghost) #variable for state's legal actions
            topScore = float("inf") #Initialize best score to +infinity
            curScore = topScore #Holds our current score, initialized equal to top score
            #Base Cases
            if gameState.isWin() or gameState.isLose():
                return gameState.getScore()
            if ghost == gameState.getNumAgents()-1:
                nxtGhost = 0
            #Begin Loop
            for action in actions:
                #If last ghost in ghosts' turn
                if nxtGhost == 0:
                    if depth == self.depth - 1:
                        curScore = self.evaluationFunction(gameState.generateSuccessor(ghost, action))
                    else:
                        curScore = minimaxDispatch(gameState.generateSuccessor(ghost, action), depth+1, alpha, beta)
                #otherwise we recursively call our max function
                else:
                    curScore = maxValue(gameState.generateSuccessor(ghost, action), depth, nxtGhost, alpha, beta)
                #if our current score is higher than our best score, replace
                if curScore < topScore:
                    topScore = curScore
                #update beta while topscore changes
                beta = min(beta, topScore)
                if topScore < alpha:
                    return topScore
            #Return current best score
            return topScore
        #Recursively call our dispatch function
        return minimaxDispatch(gameState,0, float("-inf"), float("inf"))

class ExpectimaxAgent(MultiAgentSearchAgent):
    """
      Your expectimax agent (question 4)
    """

    def getAction(self, gameState):
        """
          Returns the expectimax action using self.depth and self.evaluationFunction
          All ghosts should be modeled as choosing uniformly at random from their
          legal moves.
        """
        #Declarations
        actions = gameState.getLegalActions(0) #variable for state's legal actions
        topScore = float('-inf') #Initialize best score to -infinity
        topAction  = Directions.STOP #Initialize best action to neutral STOP
        #Begin dispatch loop
        #looping for best score from successors
        for action in actions:
            curScore = self.minValue(gameState.generateSuccessor(0, action), 1, 0)
            # if our current score is higher than our best score, replace
            if ((curScore) > topScore ):
                topScore = curScore
                #Also update top action
                topAction = action
        #Return final action taken
        return topAction
    def maxValue(self, gameState, depth):
        """For the Max Player here Pacman"""
        #Declaration
        actions = gameState.getLegalActions(0) #variable for state's legal actions
        #Begin recursion
        if ((depth == self.depth)  or (len(actions) == 0)):
            return self.evaluationFunction(gameState)
        return max([self.minValue(gameState.generateSuccessor(0, action), 1, depth) for action in actions])

    def minValue(self, gameState, ghost, depth):
        """ For the MIN Players or Agents  """
        #Declarations
        actionCount = len(gameState.getLegalActions(ghost)) #variable number of actions taken
        #Begin Recursion
        #If there are not legal actions remaining
        if (actionCount == 0):
            return self.evaluationFunction(gameState)
        #if our agent index is less than number of actions in state
        if (ghost < gameState.getNumAgents() - 1):
            return sum([self.minValue(gameState.generateSuccessor(ghost, action), ghost + 1, depth) for action in gameState.getLegalActions(ghost)]) / float(actionCount)
        #else we have reached the final ghost
        else:  #the last ghost HERE
            return sum([self.maxValue(gameState.generateSuccessor(ghost, action), depth + 1) for action in gameState.getLegalActions(ghost)]) / float(actionCount)



def betterEvaluationFunction(currentGameState):
    """
    Your extreme ghost-hunting, pellet-nabbing, food-gobbling, unstoppable
    evaluation function (question 5).

    DESCRIPTION: <write something here so we know what you did>
    """
    "*** YOUR CODE HERE ***"
    util.raiseNotDefined()

# Abbreviation
better = betterEvaluationFunction
