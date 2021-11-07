import sys
sys.path.append(".")

# Import the student solutions
import Lab7

import pathlib
DIR=pathlib.Path(__file__).parent.absolute()

import joblib
answers = joblib.load(str(DIR)+"/answers_Lab7.joblib")
print("Keys",answers.keys())

import numpy as np
import pandas as pd
answer_tol = 1e-10

def test_exercise_1():
    num_blocks = len(Lab7.lens)
    assert num_blocks == answers['exercise_1_num_blocks']
    lower = np.all(Lab7.counts > answers['exercise_1_means'] - 3*answers['exercise_1_stdevs'])
    assert lower
    upper = np.all(Lab7.counts < answers['exercise_1_means'] + 3*answers['exercise_1_stdevs'])
    assert upper

def test_exercise_2():
    answer = Lab7.greedy_sorting(Lab7.P)
    assert answers["exercise_2"] == answer

def test_exercise_3():
    P_list2 = [3,4,5,-12,-8,-7,-6,1,2,10,9,-11,13,14]
    P2 = pd.Series(P_list2,index=list(range(1,len(P_list2)+1)))
    nbreakpoints_P2 = Lab7.count_breakpoints(P2)
    P_list3 = [3,4,5,-12,-8,-7,-6,1,2,10,9,-11,14,13]
    P3 = pd.Series(P_list3,index=list(range(1,len(P_list2)+1)))
    nbreakpoints_P3 = Lab7.count_breakpoints(P3)
    assert answers["exercise_3_nbreakpoints_P2"] == nbreakpoints_P2
    assert answers["exercise_3_nbreakpoints_P3"] == nbreakpoints_P3

def fix_edge_list(edge_list):
    for i in range(len(edge_list)):
        edge_list[i] = tuple(np.sort(edge_list[i]))
    return edge_list

def test_exercise_4():
    G = Lab7.genome_to_graph([pd.Series([1,-2,-3,4]),pd.Series([5,6,7,8,9,10])])
    edge_list = fix_edge_list(Lab7.to_edge_list(G))
    answers['exercise4_edge_list'] = fix_edge_list(answers['exercise4_edge_list'])
    assert set(answers['exercise4_edge_list']) == set(edge_list)

def test_exercise_5():
    P4_list = [1,-2,-3,4]
    P4 = pd.Series(P4_list)
    P5_list = [1,3,2,-4]
    P5 = pd.Series(P5_list)

    G_P4_P5 = Lab7.combine(Lab7.genome_to_graph([P4]),Lab7.genome_to_graph([P5]))
    answers['exercise5_edge_list'] = fix_edge_list(answers['exercise5_edge_list'])
    edge_list = fix_edge_list(Lab7.to_edge_list(G_P4_P5))    
    assert set(answers['exercise5_edge_list']) == set(edge_list)

def test_exercise_6():
    P4_list = [1,-2,-3,4]
    P4 = pd.Series(P4_list)
    P5_list = [1,3,2,-4]
    P5 = pd.Series(P5_list)
    ncycles = Lab7.cycles(Lab7.genome_to_graph([P4]),Lab7.genome_to_graph([P5]))
    assert answers['exercise6_ncycles'] == ncycles

def test_exercise_7():
    P6_list = [1,2,3,4,5,6]
    P6 = pd.Series(P6_list)
    P7_list = [1,-3,-6,-5]
    P7 = pd.Series(P7_list)
    P8_list = [2,-4]
    P8 = pd.Series(P8_list)
    distance = Lab7.two_break_distance(Lab7.genome_to_graph([P6]),Lab7.genome_to_graph([P7,P8]))
    assert answers['exercise7_distance'] == distance

def test_exercise_8():
    test_edge_cycle = [[1, -3], [-3, -4], [-4, -1], [-1, 4], [4, 2], [2, 1]]
    checked_cycle, colors = Lab7.red_blue_cycle_check(Lab7.G_P4_P5,test_edge_cycle)
    assert np.all(answers['exercise8_colors'] == colors)

def test_exercise_9():
    steps = Lab7.shortest_rearrangement_scenario([pd.Series([1,-2,-3,4])],[pd.Series([1,2,-4,-3])])
    assert answers['exercise9_last_step'] == steps[-1]
