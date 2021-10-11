# Complete the sorting algorithms below.
# Since Python uses pass-by-reference for function calls
# with objects, you should sort the copy of the input
# list instead of the input list itself.

# Question 1: Bubble Sort
from math import floor


def bubbleSort(inputLst):
    lst = inputLst.copy()
    n = len(lst)
    start = 0
    fin = n-1
    hold = 0
    swapped = 1
    # Begin sort
    while swapped != 0:
        swapped = 0
        # checking each element for out of order elements
        for i in range(start, fin):
            if lst[i] > lst[i+1]:
                lst[i], lst[i+1] = lst[i+1], lst[i]
                swapped = 1
    # While loop terminates here when swapped is false
    return lst

# Question 2: Cocktail Shaker Sort
def cocktailShakerSort(inputLst):
    lst = inputLst.copy()
    n = len(lst)
    swapped = 1
    while not(swapped == 0):
        swapped = 0
        for i in range(0, n - 2):
            if lst[i] > lst[i+1]: #If elements are in wrong order
                lst[i], lst[i+1] = lst[i+1], lst[i] #Swap elements
                swapped = 1
            #end if
        #end for
        if swapped == 0:
            break
        #end if
        for i in range(n - 2, 0, -1):
            if lst[i] > lst[i+1]: #elements in wrong order
                lst[i], lst[i+1] = lst[i+1], lst[i] #swap elements
                swapped = 1
            #end if
        #end for
    #end loop
    return lst

# Question 3: Comb Sort
def combSort(inputLst):
    lst = inputLst.copy()
    n = len(lst)
    gap = n #initialize gap size
    shrink = 1.3 #Gap shrink factor
    sorted = 0
    while sorted == 0:
        #updating dap value for a next comb
        gap = floor(gap / shrink)
        if gap <= 1:
            gap = 1
            sorted = 1
        #end if
        i = 0
        while i + gap < len(lst):
            if lst[i] > lst[i+gap]:
                lst[i], lst[i+gap] = lst[i+gap], lst[i] #performing swap
                sorted = 0
            #end if
            i = i + 1
        #end loop
    #end loop
    return lst

# Question 4: Gnome Sort
def gnomeSort(inputLst):
    lst = inputLst.copy()
    n = len(lst)
    pos = 0
    while pos < n:
        if pos == 0 or lst[pos] >= lst[pos-1]:
            pos = pos + 1
        #end if
        else:
            lst[pos], lst[pos-1] = lst[pos-1], lst[pos] #performing swap
            pos = pos -1
        #end else
    #end loop
    return lst
