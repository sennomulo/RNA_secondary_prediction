# 更新日期：2022.3.18

# # 数据类型说明
# 1. sequence：list 形式，`['A', 'C', 'G', 'U',...]`
# 2. pair_list: `[[first, end, [k1, k2], [first, end, [k1, k2],...]]`
# 3. pair: `[first, end, k]`
# 4. state: `[[first, end, k], [first, end, k],...]`

import numpy as np
# import pandas as pd


def single_pair(sequence, a,b):
    condition1 = (sequence[a] == 'U' and sequence[b]=='A') or (sequence[a]=='A' and sequence[b]=='U') 
    condition2 = (sequence[a]=='C' and sequence[b]=='G') or (sequence[a]=='G' and sequence[b]=='C')
    condition3 = (sequence[a] == 'G' and sequence[b] == 'U') or (sequence[a]=='U' and sequence[b]=='G')
    return (condition1 or condition2 or condition3)



def consecutive_base_pairs(sequence, min_stem, min_loop):
    # 这种办法应该算是有效率的了，
    # 如果是按 minStem 和 minLoop 的长度「批量」检测的话，
    # 就不好搜索到大于 minStem 和 minLoop 的碱基对。
    pairs = []
    n = len(sequence)
    for i in range(n-2*min_stem - min_loop+1):
        for j in range(i+2*min_stem + min_loop -1, n):
            con_pair = 0
            k = 0
            while (k<= (j-i-min_loop+1)/2):
                if single_pair(sequence, i+k, j-k):
                    con_pair+=1
                    if con_pair >= min_stem:
                        temp_pair = [i,j,[con_pair]]
                        pairs.append(temp_pair)
                else:
                    break
                
                k+=1
    return pairs


def K_consecutive_base_pairs(pairs):
    pair_list = []
    temp_rd1 = ' '
    temp_rd2 = ' '

    for i in range(len(pairs)):
        if (pairs[i][0] == temp_rd1 and pairs[i][1]==temp_rd2):
            pair_list[len(pair_list)-1][2].append(pairs[i][2][0])
        else:
            pair_list.append(pairs[i])
            temp_rd1 = pairs[i][0]
            temp_rd2 = pairs[i][1]

    return pair_list


def random_choose(pair_list):
    '''
    用来随机挑选 pair_list 里的一个碱基对\n
    返回一个 [i, j, k]
    '''
    i = np.random.randint(len(pair_list))
    k = np.random.randint(len(pair_list[i][2]))
    chosen_pair = [pair_list[i][0], pair_list[i][1], pair_list[i][2][k]]
    return chosen_pair

def specific_pairs(pair):
    '''
    给出 pair:[a,b,c] 的具体碱基配对情况（只告诉序号）\n
    返回一个[[a, b], [a, b], ......] 和 [一堆数字]
    '''
    pairs_coupled = []
    pairs_order = []
   
   
    for k in range(pair[2]):
        pairs_order.append(pair[0]+k)
        pairs_order.append(pair[1]-k)
        pairs_coupled.append([pair[0]+k, pair[1]-k])
      

    return pairs_coupled, pairs_order


def remove_conflict(initial_state, taboo, neighbor_state, min_stem)->None:
    '''
    输入初态、冲突集和邻态，就可以把初态和冲突集的冲突部分删去，
    再把符合要求的补到邻态中。
    '''
    for pair in initial_state:
        remain = []
        pairs_couple, _ = specific_pairs(pair)
        for subpair in pairs_couple:
            check = True
            for taboo_num in taboo:
                if taboo_num in subpair:
                    check = False
            if check:
                remain.append(subpair)

                
        
        if len(remain) >= min_stem:
            start = remain[0][0]
            end = remain[0][1]
            length = len(remain)
            neighbor_state.append([start, end, length])

def find_neighbor(pair_list, initial_state, min_stem):
    '''
    寻找一个邻态\n
    返回一个[[i, j, k], [i, j, k], ...]
    '''
    second_pair = random_choose(pair_list)
    neighbor_state = []
    neighbor_state.append(second_pair)

    # 冲突集
    _, taboo = specific_pairs(second_pair)
    
    # 挨个查，删冲突
    remove_conflict(initial_state, taboo, neighbor_state, min_stem)
    
    return neighbor_state

def total_pairs(state):
    '''
    计算序列中所有碱基对的数量
    '''
    TP = 0
    for pair in state:
        TP += pair[2]
    return TP


def pseudoknots(state):
    '''
    计算假结数量
    '''
    count = 0
    for i in range(len(state)-1):
        for j in range(i+1, len(state)):
            cond1 = state[i][0] < state[j][0] < state[i][1] < state[j][1]
            cond2 = state[j][0] < state[i][0] < state[j][1] < state[i][1]
            if cond1 or cond2:
                count += 1
    return count


def fitness(state, MG)->float:
    """
    计算适应度函数
    ---
    输入：\n
    state: 状态[[a,b,c], [a,b,c], ...]\n
    MG: 期望的假结数量
    """
    Group = len(state)
    TP = total_pairs(state)
    AP = float(TP)/float(Group)
    PG = pseudoknots(state)

    if PG <= MG:
        return TP * AP **2
    else:
        return TP * AP **2 * (Group-PG) / Group


def energy_data(neighbor_state, max_state, MG):
    max_fitness = fitness(max_state, MG)
    neighbor_fitness = fitness(neighbor_state, MG)
    return (neighbor_fitness-max_fitness)


def framework(sequence, min_stem, min_loop, max_T, min_T, T_step, MG, epoch):
    '''
    输出一个 max_state

    输入值如下：
    sequence: 序列
    min_stem: 最小茎长度
    min_loop: 最小环长度
    max_T: 最高温度
    min_T: 最低温度
    T_step: 温度下降值
    MG: 预期的假结数量
    epoch: 每一个温度下的循环次数
    '''
    # 初始工作
    max_state = []
    T = max_T

    # 寻找连续碱基对集
    pairs = consecutive_base_pairs(sequence, min_stem, min_loop)
    pair_list = K_consecutive_base_pairs(pairs)

    # 随机挑选初态
    max_state.append(random_choose(pair_list))

    while(T>min_T):
        for _ in range(epoch):
            neighbor_state = find_neighbor(pair_list, max_state, min_stem)
            dE = energy_data(neighbor_state, max_state, MG)
            if dE >= 0:
                max_state = neighbor_state
            else:
                if np.exp(dE/T) > np.random.random():
                    max_state = neighbor_state
        T -= T_step * np.sqrt(T) # 让高温降得更快，低温降得更慢

    return max_state


def recognize_symbols(first, pair_position):
    '''
    输入：
    first: 第一个碱基的序列
    pair_position: [[left first, left end, right end, right first], [], ...]\n
    说白了就是根据网站上的格式来进行的
    
    输出：所包含的所有 pair 和假结数量。
    '''
    state = []
    pairs_coupled = []
    for pair in pair_position:
        state.append([pair[0]-first, pair[3]-first, pair[1]-pair[0]+1])
        for k in range(pair[1]-pair[0]+1):
            pairs_coupled.append([pair[0]-first+k, pair[3]-first-k])
    
    MG = pseudoknots(state)
    return pairs_coupled, MG


def to_base_pairs(state):
    '''
    给出 state，返回 state 中所包含的所有 pair.
    '''
    pairs_coupled = []
    for pair in state:  
        for k in range(pair[2]):           
            pairs_coupled.append([pair[0]+k, pair[1]-k])
      
    return pairs_coupled


def evaluation(actual, prediction):
    '''
    输入真实的和预测的碱基对，输出 SP、SN 和 F_measure (单位均为 %)
    '''
    prediction_length = len(prediction)
    actual_length = len(actual)
    
    TP = 0
    for pair in prediction:
        if pair in actual:
            TP+=1
    FP = prediction_length - TP
    FN = actual_length - TP

    SN = TP/ (TP+FN) *100
    SP = TP/(TP+FP) *100
    F_measure = 2*SN*SP/(SN+SP+0.0001) # 小数用于防止分母为零
    return SP,SN,F_measure


# 全自动估计
def auto_eval(sequence, pairing, min_stem=2, min_loop=3, max_T=100, min_T=20, T_step=5, epoch=20):
    '''
    返回 max_state_interpretation, (SP,SN,F_measure)
    '''
    
    # pairing = [first, [], [], ...]
    actual, MG = recognize_symbols(pairing[0], pairing[1:])
    max_state = framework(sequence, min_stem, min_loop, max_T, min_T, T_step, MG, epoch)
    prediction = to_base_pairs(max_state)
    SP,SN,F_measure = evaluation(actual, prediction)

    # 转化 max_state
    max_state_interpretation = interpretation(prediction,len(sequence),pairing[0])
    return max_state_interpretation, (SP,SN,F_measure)
    
def interpretation(pairs_coupled, length, first):
    original = np.arange(length)
    for pair in pairs_coupled:
        original[pair[0]] = pair[1]
        original[pair[1]] = pair[0]
    interpret = original + first
    return interpret.tolist()