import math
def diff_o_sum(test_A, test_B)
    diff=sum([a_i - b_i for a_i, b_i in zip(test_A, test_B)])
    return diff
    
def labe_similarity(v_A,v_B):#CREATE VECTORS USING THE SELECTED TERMS...
    "compute labe similarity"
    test_A=v_A
    test_B=v_B
    len_A=len(test_A)
    len_B=len(test_B)
