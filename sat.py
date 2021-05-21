import sympy as sp

def Sat(agent, prop):
	return [i for i in range(n_states) if prop in label[agent][i]]

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
	return (len(set(a).union(set(b))) - len(set(a).intersection(set(b))))/len(set(a).union(set(b)))

def sat(epsilon):
	sat = []
	PhiM = {}
	Phie = {}
	for s in states:
		PhiM[s] = []
		Phie[s] = []
		for s1 in Sat("M","p"):
			if h("M",[s],s1) == 1:
				PhiM[s].append(s1)
		for s1 in Sat("e","p"):
			if h("e",[s],s1) == 1:
				Phie[s].append(s1)

	for s in states:
		if epsilon["accuracy"](1 - jaccard(PhiM[s],Phie[s])):
			sat.append(s)

	return sat

if __name__=="__main__":
	# Model
	global states
	global agents
	global transition
	global probability
	global ap
	global label
	global n_states

	n_states = 16
	states = range(n_states)
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

	# CHIEDERE AD ALBERTO
	probability = {}
	probability["M"] = [0, 1, 0, 0]
	probability["e"] = [0, 0.5, 0.5, 0]

	ap = {}
	ap["M"] = ["p","q","r","s"]
	ap["e"] = ["p","q","r","s"]

	label = {}

	label["M"] = {}
	label["M"][0] = {}
	label["M"][1] = {"p"}
	label["M"][2] = {"q"}
	label["M"][3] = {"r"}
	label["M"][4] = {"s"}
	label["M"][5] = {"p","q"}
	label["M"][6] = {"p","r"}
	label["M"][7] = {"p","s"}
	label["M"][8] = {"q","r"}
	label["M"][9] = {"q","s"}
	label["M"][10] = {"r","s"}
	label["M"][11] = {"p","q","r"}
	label["M"][12] = {"p","q","s"}
	label["M"][13] = {"p","r","s"}
	label["M"][14] = {"q","r","s"}
	label["M"][15] = {"p","q","r","s"}

	label["e"] = {}
	label["e"][0] = {}
	label["e"][1] = {"p"}
	label["e"][2] = {"q"}
	label["e"][3] = {"r"}
	label["e"][4] = {"s"}
	label["e"][5] = {"p","q"}
	label["e"][6] = {"p","r"}
	label["e"][7] = {"p","s"}
	label["e"][8] = {"q","r"}
	label["e"][9] = {"q","s"}
	label["e"][10] = {"r","s"}
	label["e"][11] = {"p","q","r"}
	label["e"][12] = {"p","q","s"}
	label["e"][13] = {"p","r","s"}
	label["e"][14] = {"q","r","s"}
	label["e"][15] = {"p","q","r","s"}

	# Formula
	epsilon = {}
	# epsilon["formula"] = {}
	# epsilon["formula"]["boolean"] = [True, "p"]
	# epsilon["formula"]["temporal"] = lambda x,y:
	# epsilon["formula"][""] = 6
	epsilon["probability"] = lambda x: x <= 0.5
	epsilon["accuracy"] = lambda x: x <= 0.6

	print(sat(epsilon))