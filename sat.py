import sympy as sp

def sat_agent(agent, phi, until=False, next=False):

    if next:
        custom_hitting_function = lambda t_a, t_b, s : h_next(agent,t_a,s)
    elif (not until):
        custom_hitting_function = lambda t_a, t_b, s : h_cond(agent,t_a,t_b,s)
    else:
        custom_hitting_function = lambda t_a, t_b, s : h_cond_bounded(agent,t_a,t_b,s,until)

    temp_before = set()

    for disjunct in phi["psi"]["before"]:
        temp = set(states)
        for conjunct in disjunct:
            if conjunct[0] != "-":
                temp = temp.intersection([i for i in states if conjunct in label[agent][i]])
            else:
                temp = temp.intersection([i for i in states if conjunct[1:] not in label[agent][i]])
        temp_before = temp_before.union(temp)

    temp_after = set()

    for disjunct in phi["psi"]["after"]:
        temp = set(states)
        for conjunct in disjunct:
            if conjunct[0] != "-":
                temp = temp.intersection([i for i in states if conjunct in label[agent][i]])
            else:
                temp = temp.intersection([i for i in states if conjunct[1:] not in label[agent][i]])
        temp_after = temp_after.union(temp)

    return set([ s for s in states if phi["probability"]( custom_hitting_function(temp_after,temp_before,s) ) ])

def h_cond(agent,A,B,initial_state):
    syms = sp.symbols("x0:"+str(n_states))
    variables = sp.Matrix(syms[0:n_states])
    setOfEquations = sp.Matrix(transition[agent]) * variables

    for i in A:
        setOfEquations.row_del(i)
        setOfEquations = setOfEquations.row_insert(i, sp.Matrix([1]))
        
    for i in (states-A-B):
        setOfEquations.row_del(i)
        setOfEquations = setOfEquations.row_insert(i, sp.Matrix([0]))

    eq = sp.Eq(setOfEquations,variables)
    ans = (sp.solve(eq,syms,dict=True))[0].get(syms[initial_state])

    if ans:
        return ans.subs((x,0) for x in syms)
    else:
        return 0
    
def h_cond_bounded(agent,A,B,initial_state,t):
    modified_transition_matrix = sp.Matrix(transition[agent])
    
    for i in states-B:
        modified_transition_matrix.row_del(i)
        modified_transition_matrix = modified_transition_matrix.row_insert(i, sp.zeros(1,n_states))
        modified_transition_matrix[i,i] = 1
    
    modified_transition_matrix = modified_transition_matrix**t

    sum = 0
    for i in A:
        sum = sum + modified_transition_matrix[initial_state,i]

    return sum

def h_next(agent,A,initial_state):
    transition_matrix = sp.Matrix(transition[agent])

    sum = 0
    for i in A:
        sum = sum + transition_matrix[initial_state,i]

    return sum
    
def h(agent,A,initial_state):
    syms = sp.symbols("x0:"+str(n_states))
    variables = sp.Matrix(syms[0:n_states])
    setOfEquations = sp.Matrix(transition[agent]) * variables

    for i in A:
        setOfEquations.row_del(i)
        setOfEquations = setOfEquations.row_insert(i, sp.Matrix([1]))

    eq = sp.Eq(setOfEquations,variables)
    ans = (sp.solve(eq,syms,dict=True))[0].get(syms[initial_state])

    if ans:
        return ans.subs((x,0) for x in syms)
    else:
        return 0

def jaccard(a,b):
    if ( set(a).union(set(b)) == set()):
        return 0
    else:
        return 1 - (len(set(a).union(set(b))) - len(set(a).intersection(set(b))))/len(set(a).union(set(b)))

def sat(epsilon):
    
    out = []
    PhiM = {}
    Phie = {}
    for s in states:
        PhiM[s] = []
        Phie[s] = []

        for s1 in sat_agent("M",epsilon["phi"]):
            if h("M",[s1],s) == 1:
                PhiM[s].append(s1)

        for s1 in sat_agent("e",epsilon["phi"]):
            if h("e",[s1],s) == 1:
                Phie[s].append(s1)

    for s in states:
        ji = jaccard(PhiM[s],Phie[s])
        print("State "+str(s)+" has jaccard index = "+str(ji))
        if epsilon["accuracy"](ji):
            out.append(s)

    return set(out)

if __name__=="__main__":
    global states
    global agents
    global transition
    global probability
    global ap
    global label
    global n_states
    
    # Representation of the models "M" and "e" of the example in
    # "Modelling accuracy and trustworthiness of explaining agents"
    # submitted to the conference LORI 2021

    n_states = 16
    states = set(range(n_states))
    agents = {"M","e"}

    transition = {}
    transition["M"] = [[0.25, 0.25, 0.25, 0.25, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 1   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0.25, 0   , 0   , 0.25, 0.25, 0.25, 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 ]]

    transition["e"] = [[0.25, 0.25, 0.25, 0.25, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0.5 , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0.25, 0   , 0   , 0.25, 0.25, 0.25, 0   , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.25, 0.25, 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   , 0   , 0   , 0.5 , 0   , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 , 0   ],
                      [0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0.5 , 0.5 ]]

    label = {}

    label["M"] = {}
    label["M"][0] = {}
    label["M"][1] = {"r"}
    label["M"][2] = {"q"}
    label["M"][3] = {"r"}
    label["M"][4] = {"s"}
    label["M"][5] = {"r","q"}
    label["M"][6] = {"q","r"}
    label["M"][7] = {"q","s"}
    label["M"][8] = {"q","r"}
    label["M"][9] = {"q","s"}
    label["M"][10] = {"r","s"}
    label["M"][11] = {"q","r"}
    label["M"][12] = {"p","q","s"}
    label["M"][13] = {"r","s"}
    label["M"][14] = {"q","r","s"}
    label["M"][15] = {"q","r","s"}

    label["e"] = {}
    label["e"][0] = {}
    label["e"][1] = {"r"}
    label["e"][2] = {"q"}
    label["e"][3] = {"r"}
    label["e"][4] = {"s"}
    label["e"][5] = {"r","q"}
    label["e"][6] = {"q","r"}
    label["e"][7] = {"q","s"}
    label["e"][8] = {"q","r"}
    label["e"][9] = {"q","s"}
    label["e"][10] = {"r","s"}
    label["e"][11] = {"q","r"}
    label["e"][12] = {"p","q","s"}
    label["e"][13] = {"r","s"}
    label["e"][14] = {"q","r","s"}
    label["e"][15] = {"q","r","s"}

    # Formula
    # A^{e}_{>= 0.75} P_{=1} (\top \cup p)
    epsilon = {}
    epsilon["phi"] = {}
    epsilon["phi"]["probability"] = lambda x: x == 1
    epsilon["phi"]["psi"] = {}
    epsilon["phi"]["psi"]["before"] = [[ ]] # [[ ]] stands for \top
    epsilon["phi"]["psi"]["after"] = [[ "p" ]]
    epsilon["accuracy"] = lambda x: x >= 0.75
    
    # TODO: Extend epsilon to most general form:
    #
    # epsilon := *A^{Gamma}_{nabla h} phi
    # phi := p |
    #        phi1 and phi2 |      // sat(phi1) intersect sat(phi2)
    #        neg phi |            // complement of sat(phi)
    #        *P_{nabla p} psi
    #
    # psi := *next phi |
    #        *phi until phi |
    #        *phi until_{<= t} phi
    #
    # * = done

    print(sat_agent("M",epsilon["phi"],next=True))

    s = sat(epsilon)
    
    print("*****")
    print("The following states satisfy the accuracy requirement: " + str(s))
    