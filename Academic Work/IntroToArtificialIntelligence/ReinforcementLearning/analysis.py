# analysis.py
# -----------
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


######################
# ANALYSIS QUESTIONS #
######################

# Set the given parameters to obtain the specified policies through
# value iteration.

def question2():
    answerDiscount = 0.9
    answerNoise = 0.01
    #Decreasing our noise allows the agent to more successfully make intended actions
    #and thus cross the bridge more successfully
    return answerDiscount, answerNoise

#prefer close exit (+1), risking the cliff (-10)
def question3a():
    #discount of .1 because since we are risking the cliff it should be getting a lot smaller
    answerDiscount = .1
    #no noise because we are risking the cliff so there is a 0% chance this goes as planned
    answerNoise = 0
    #no living reward because we are risking the cliff and going for the close exit
    answerLivingReward = 0
    return answerDiscount, answerNoise, answerLivingReward
    # If not possible, return 'NOT POSSIBLE'
#prefer close exit (+1), but avoiding the cliff (-10)
def question3b():
    #discount of .2 because we are avoid the cliff but still preferring the close exit
    #so it should be pretty small but not as small as risking the cliff
    answerDiscount = .2
    #.2 noise because there is a 20% chance this goes as planned since we are avoiding the cliff
    answerNoise = .2
    #negative reward because we avoid the cliff but are only going for the close exit
    answerLivingReward = -1
    return answerDiscount, answerNoise, answerLivingReward
    # If not possible, return 'NOT POSSIBLE'
#prefer distant exit (+10) risking the cliff (-10)
def question3c():
    #discount of .9 because we want the higher reward exit so it will
    answerDiscount = .9
    #no noise because we are risking the cliff so there is a 0% chance this goes as planned
    answerNoise = 0
    #negative reward because we have to go all the way to the distant exit while risking the cliff
    answerLivingReward = -1
    return answerDiscount, answerNoise, answerLivingReward
    # If not possible, return 'NOT POSSIBLE'
#prefer distant exit (+10), avoiding cliff (-10)
def question3d():
    #higher discount than 3c because we are avoiding the cliff rather than risking it
    answerDiscount = 1
    #.4 noise because there is a 40% chance this goes as planned since we are avoiding the cliff
    #but also going for thte distant exit
    answerNoise = .4
    #no reward because we completely avoid the cliff but are also going for the far exit
    answerLivingReward = 0
    return answerDiscount, answerNoise, answerLivingReward
    # If not possible, return 'NOT POSSIBLE'
#avoid both exits AND cliff
def question3e():
    #no discount because we are avoiding everything
    answerDiscount = 0
    #no noise because we are risking the cliff so there is a 0% chance this goes as planned
    answerNoise = 0
    #no reward because avoiding everything
    answerLivingReward = 0
    return answerDiscount, answerNoise, answerLivingReward
    # If not possible, return 'NOT POSSIBLE'

def question8():
    answerEpsilon = None
    answerLearningRate = None
    return answerEpsilon, answerLearningRate
    # If not possible, return 'NOT POSSIBLE'

if __name__ == '__main__':
    print('Answers to analysis questions:')
    import analysis
    for q in [q for q in dir(analysis) if q.startswith('question')]:
        response = getattr(analysis, q)()
        print('  Question %s:\t%s' % (q, str(response)))
