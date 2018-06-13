import random
import copy

print gap('LoadPackage("numericalsgps");')

def ConvertGapToSage(var):
	try:
		ret = [ConvertGapToSage(j) for j in var]

		return ret
	except:
		pass

	# bool conversion
	if var == gap('true') or var == gap('false'):
		ret = bool(var)

		return ret

	# int conversion
	try:
		ret = Rational(var)

		if ret.denominator() == 1 :
			ret = int(ret)

		return ret
	except:
		pass

	return 0 


# Helper Functions -- Plotting **********************************************************************************************

# general plotting function
#	input:	list of points, values for each point
#	output: 2d plot
def GeneralPlotting(pts, vals):
	if len(pts) == len(vals):
		plot = Graphics()
		plot += list_plot(pts, size = 25, color = 'grey')
		for i in range(len(pts)):
			plot += text('     '+ str(vals[i]), pts[i], fontsize = 'small', color = 'blue', figsize = (15,15)) 
		
		return plot

# highlights pts on plot 
def CirclePoints(plot, pts, color):
	ret = Graphics()
	ret += plot
	ret += list_plot(pts, alpha = 150, size = 150, color = color)

	return ret


# Helper Functions -- Semigroup **********************************************************************************************

# creates affine semigroup in GAP
def CreateAffineSemigroup(gens,name=""):
	genstr = ','.join(str(x) for x in sorted(gens))
	if name == "":
		return gap('AffineSemigroupByGenerators(' + genstr + ');')
	else:
		gap.eval(name + ' := AffineSemigroupByGenerators(' + genstr + ');')
		return gap(name + ';')

# prints n random generators in the range of (1,10) 
def GenerateGenerators(nmax,ngens):
	return list([random.randint(1,nmax),random.randint(1,nmax)] for n in range(ngens))

# generating points in semigroup from (s,s) to (n,n):
def GeneratePointsFrom(semigroup,n,s):
	return [[x,y] for x in range(s,n+1) for y in range(s,n+1) if [x,y] in semigroup]

# generating max factorization length for the pts given to the corresponding given semigroup:
def GenerateMaxFact(semigroup,pts):
	semigroup.FactorizationsUpToElement(max(pts))
	return [max(semigroup.LengthSet(pt)) for pt in pts]

# generating number of factorizations for the pts given to the corresponding given semigroup:
def GenerateNumFact(semigroup,pts):
	semigroup.FactorizationsUpToElement(max(pts))
	return [len(semigroup.Factorizations(pt)) for pt in pts]

# finding the outer generators of 3 generator semigroup
def FindOuterGenerators(semigroup):
	minPres = semigroup.MinimalPresentation()
	for p in minPres[0]:
		if (p[0] == 0 and p[1] == 0) or (p[1] == 0 and p[2] == 0) or (p[0] == 0 and p[2] == 0):
			return p

	return [0,0,0]


# Helper Functions -- Decomposition **********************************************************************************************

# number of factorizations
def NumFact(group, pt):
	return len(group.Factorizations(pt))

# max number of factorizations
def MaxFact(group, pt):
	return max(group.LengthSet(pt))

# min number of factorizations
def MinFact(group, pt):
	return min(group.LengthSet(pt))

def RemainingValPts(pts, vals):
	ret_pts = []
	for pt in pts:
		if vals[pts.index(pt)] != 0:
			ret_pts.append(pt)

	return ret_pts

def SuperImposeConeVals(coneList, pts, vals):
	ret_pts = list(pts)
	ret_vals = list(vals)
	for pt in pts:
		i = pts.index(pt)
		for cone in coneList:
			if cone.ContainsPoint(pt):
				ret_vals[i] = ret_vals[i] - cone.func(cone.group, pt)#cone.NumFact(pt)
	return ret_vals

def isInt(input):
    try: int(input)
    except ValueError: return False
    else: return True


# Cone Decomposition Class ******************************************************************************************************

class ConeDecomposition:
	def __init__(self, semigroup, xmax, ymax, func=None, degree=None): 
		self.semigroup = semigroup
		self.group_gens = semigroup.gens
		self.MinimalPresentation = semigroup.MinimalPresentation()
		self.ElementofMinPres = semigroup.ElementsOfMinimalPresentation()
		self.semigroup.FactorizationsUpToElement([xmax,ymax])

		self.degree = degree

		if self.degree == None:
			self.degree = 1

		if func == None:
			func = lambda x,y: 1
		
		self.func = func

		self.label_color = ["red","orange","green","blue","purple","aqua","magenta","olive","black","brown","coral","crimson","cyan","gold","greenyellow","mistyrose","salmon"]

		self.__decompositions = {}
		#self.Decomposition(xmax,ymax)

	def __InitConeList(self):
		self.nCones = 0
		self.coneList = []
		self.plot = copy.copy(self.group_plot)

	def FindDecomposition(self,xmax,ymax):
		#if self.__decompositions[str([xmax,ymax])] != None:
		if str([xmax,ymax]) in self.__decompositions: 
			return self.__decompositions[str([xmax,ymax])]
		else:
			print "Decomposition not found"

	def Decomposition(self, xmax,ymax):

		# members of the affine semigroup
		self.pts = GeneratePoints(self.semigroup,xmax,ymax)
		# vals of members corresponding to func
		self.vals = [self.func(self.semigroup,pt) for pt in self.pts]
		# plot of semigroup 
		self.group_plot = self.semigroup.Plot(xmax,ymax)

		self.__InitConeList()

		print "Beginning Decomposition for Affine Semigroup generate by " + str(self.group_gens)
		print "Minimal Presentation = " + str(self.MinimalPresentation) + " @ " + str(self.ElementofMinPres)

		remain = copy.copy(self.pts)

		while(len(remain)!= 0):
			print "\nDefining Cone " + str(self.nCones+1) + ":"


			## handling generators ##
			print "	Cone generators Options- ",
			for i in range(len(self.group_gens)):
				print str(i) + ":" + str(self.group_gens[i]) + " ",
			print " "


			# takes in generators in a list form where integers indicate the generators previous defined in group or enter unique vectors
			curConeGen = input("	Please enter generators:")
			nCurConeGen = len(curConeGen)
			# updates generators for integer values -- 
			for gen in curConeGen:
				if isInt(str(gen)):
					curConeGen[curConeGen.index(gen)] = self.group_gens[gen]
			
			## handling basepoints ##
			basepoint = remain[0]
			setBasepoint = input("	Set Basepoint: (enter vector or 1 to keep set @ " + str(basepoint) + ")  ")

			if isInt(str(setBasepoint)) == False:
				while setBasepoint not in self.semigroup:
			 		setBasepoint = input("	Error, please enter an appropriate basepoint: ")
				
				basepoint = setBasepoint
			

			## cone declaration + preview ##
			print "	Preview -- Cone generate by "+ str(curConeGen) + " with basepoint @ "+ str(basepoint)

			cone = Cone(basepoint, curConeGen, self.semigroup, self.degree, self.func)

			test = cone.PlotCone(self.plot,xmax,ymax,self.label_color[self.nCones])
			test_vals = SuperImposeConeVals(self.coneList+[cone], copy.copy(self.pts), copy.copy(self.vals))

			warning_list = []
			for i in range(len(test_vals)):
				if test_vals[i] > 0:
					test += text('      '+ str(test_vals[i]), self.pts[i], fontsize = 'small', color = 'blue', figsize = (15,15))
				elif test_vals[i] < 0:
					warning_list.append(self.pts[i])
					test += text('      '+ str(test_vals[i]), self.pts[i], fontsize = 'small', color = 'red', figsize = (15,15))

			# cone preview
			show(test)


			## print warning list of overlapping cones ##
			if len(warning_list) != 0:
				print "	*** Warning -- Cone Overlap @ " + str(warning_list)


			## handle cases -- Keep Cone, Remove Current Cone, Remove Previous Cone, End Program ##
			print "	0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program"
			options = input("	Choose option: ")
			while isInt(str(options)) == False or options < 0 or options > 4:
				options = input("	Error, please enter an appropriate option: ")

			if options == 1:
				continue
			elif options == 2:
				redoCone = input("	Redefine cone(s) starting @ cone: ")
				self.nCones = redoCone - 1
				continue
			elif options == 3:
				self.coneList.append(cone)
				print"	Cone generate by "+ str(cone.gens) + " with basepoint @ "+ str(cone.basepoint)
				print "Incomplete Decomposition finished with "+ str(self.nCones) + " cones"
				return self.coneList
				self.__decompositions[str([xmax,ymax])] = copy.copy(self.coneList)
			elif options == 4:
				return self.coneList

			self.coneList.append(cone)
			self.plot = cone.HighlightCone(self.plot,xmax,ymax, color = self.label_color[self.nCones])

			self.nCones += 1
			print"	Cone generate by "+ str(cone.gens) + " with basepoint @ "+ str(cone.basepoint)

			remain = RemainingValPts(self.pts, SuperImposeConeVals(self.coneList, copy.copy(self.pts), copy.copy(self.vals)))
			if len(remain) != 0:
				basepoint = remain[0]


		print "Finished with "+ str(self.nCones) + " cones"

		# save decomposition --
		self.__decompositions[str([xmax,ymax])] = copy.copy(self.coneList)
		return self.coneList



# Cone Class **********************************************************************************************

class Cone:
	def __init__(self, basepoint, gens, semigroup, degree=None, func=None):
		self.basepoint = basepoint
		self.gens = gens
		self.group = semigroup

		if degree == None:
			self.degree = 1
		else:
			self.degree = degree

		self.ConeMembers()
		self.initDegree()

		if func == None:
			func = lambda x,y: 1

		self.func = func
		
		self.func_coeff = self.matrix.solve_right(vector([func(semigroup, pt[0:2]) for pt in self.pts]))
		self.polynomial()

	def initDegree(self):
		ans = []

		for x in self.pts:
			cur = []
			for d in range(self.degree,-1,-1):
				for i in range(d,-1,-1):
					cur.append((x[0]**i)*(x[1]**(d-i)))
			ans.append(cur)

		self.matrix = Matrix(ans)

	def polynomial(self):
		polyx = []
		polyy = []

		for i in range(self.degree,-1,-1):
			xdegree = i 
			if xdegree > 1:
				polyx.append("x**"+str(xdegree))
			else:
				if xdegree == 1:
					polyx.append("x")
				else:
					polyx.append("")

			ydegree = self.degree-i
			if ydegree > 1:
				polyy.append("y**"+str(ydegree))
			else:
				if ydegree == 1:
					polyy.append("y")
				else:
					polyy.append("")


		polyx.append("")
		polyy.append("")

		ret = ""

		for i in range(len(self.func_coeff)):
			if i == 0:
				ret += str(self.func_coeff[0]) + polyx[0] + polyy[0]
			else:
				if self.func_coeff[i] >= 0:
					ret += " + " + str(self.func_coeff[i]) + polyx[i] + polyy[i]
				else:
					ret += " " + str(self.func_coeff[i]) + polyx[i] + polyy[i]

		self.polynomial = ret

	def ConeMembers(self):
		members = []
		curDegree = []

		for i in range(self.degree+1):
			if i == 0:
				members.append(self.basepoint)
				curDegree.append(self.basepoint)
			else:
				for d in curDegree:
					temp = []
					for g in self.gens:
						element = [d[0]+g[0],d[1]+g[1]]

						if element not in members:
							members.append(element)
							temp.append(element)

				curDegree = temp

		self.pts = members

	def DrawCone(self,plot,x_max,y_max,color):
		ret = plot 
		for gen in self.gens:
			
			if gen[0] == 0:
				x_factor = 0
			else:
				x_factor = int(x_max/gen[0])

			if gen[1] == 0:
				y_factor = 0
			else:
				y_factor = int(y_max/gen[1])

			
			max_mult = max(x_factor,y_factor)

			x = max_mult*gen[0] + gen[0] + self.basepoint[0]
			y = max_mult*gen[1] + gen[1] + self.basepoint[1]

			ret += line([self.basepoint,(x,y)],xmax = x_max, ymax = y_max, color = color)

		return ret

	def HighlightCone(self,plot,xmax,ymax,color):
		pts = [(x,y) for x in range(xmax+1) for y in range(ymax+1) if self.ContainsPoint((x,y))]
		return CirclePoints(plot, pts, color) 

	def HighlightConeWithVals(self,plot,xmax,ymax,color):
		ret = plot
		pts = [[x,y] for x in range(xmax+1) for y in range(ymax+1) if self.ContainsPoint((x,y))]
		for i in range(len(pts)):
			ret += text('        '+ str(self.func(self.group,pts[i])), pts[i], fontsize = 'small', color = 'blue', figsize = (15,15)) 
		return CirclePoints(ret, pts, color) 

	def PlotCone(self,plot,xmax,ymax,color):
		ret = self.DrawCone(plot,xmax,ymax,color)
		return self.HighlightCone(ret,xmax,ymax,color)

	def PlotConeWithVal(self,plot,xmax,ymax,color):
		ret = self.DrawCone(plot,xmax,ymax,color)
		return self.HighlightConeWithVals(ret,xmax,ymax,color)

	# only handles cases in 2 dimension
	def ContainsPoint(self,point):
		matrix = Matrix(self.gens)
		check = matrix.solve_left(vector([point[0]-self.basepoint[0],point[1]-self.basepoint[1]]))
		if check[0] >= 0 and check[1] >= 0 and check[0].denominator() == 1 and check[1].denominator() == 1:
			return True
		else:
			return False

	def __repr__(self):
		return '<Cone with basepoint %s generated by %s with polynomial function %s>' % (self.basepoint, self.gens, self.polynomial)



# Affine Semigroup Class **********************************************************************************************

class AffineSemigroup:
	def __init__(self, gens):
		if gens:
			self.__InitWithGapSemigroup(CreateAffineSemigroup(gens))

	def __InitWithGapSemigroup(self, semigroup):
		self.semigroup = semigroup
		self.gens = ConvertGapToSage(gap('GeneratorsOfAffineSemigroup('+ self.semigroup.name() + ');'))

		self.__InitWithDefaults()

		return self

	def __InitWithDefaults(self):
		self.__factorizations = {}
		self.__max = {}

	def __contains__(self,other):
		if str(other) in self.__factorizations:
			if self.__factorizations[str(other)] == []:
				return False
			else:
				return True
		else:
			if not self.BelongsToAffineSemigroup(other):
				self.__factorizations[str(other)] = []
				return False
			else:
				return True

	def GenericGapCallGlobal(self, gapfuncname):
		return ConvertGapToSage(gap(gapfuncname + '(' + self.semigroup.name() + ')'))

	# generating points in semigroup up to (xmax,ymax):
	def ElementsUpToElement(xmax,ymax):
		return [[x,y] for x in range(xmax+1) for y in range(ymax+1) if [x,y] in self]

	def Plot(self,xmax,ymax,ptsize=None,func=None):
		pts = GeneratePoints(self,xmax,ymax)
		plot = Graphics()

		if ptsize == None:
			ptsize = 25

		plot += list_plot(pts, size = ptsize, color = 'grey', figsize = (15,15))

		if func != None:
			for i in range(len(pts)):
				plot += text('        '+ str(func(self,pts[i])), pts[i], fontsize = 'small', color = 'blue', figsize = (15,15)) 

		self.plot = plot
		
		return plot
		
	# for [x,y] in Z^2
	def FactorizationsUpToElement(self, v):
		pts = sorted([[x,y] for x in range(v[0]+1) for y in range(v[1]+1)])
		self.__factorizations[str(pts[0])] = [[0 for g in self.gens]]
		self.__max[str(pts[0])] = 0

		for n in range(len(pts)):
			if pts[n] in self.semigroup:
				if str(pts[n]) in self.__factorizations:
					continue

				self.__factorizations[str(pts[n])] = []
				
				for i in range(len(self.gens)):
					temp = [x - y for (x, y) in zip(pts[n], self.gens[i])]
					if (min(temp) < 0) or (temp not in self.semigroup):
						continue

					for f in self.__factorizations[str(temp)]:
						add = list(f)
						add[i] = add[i] + 1
						self.__factorizations[str(pts[n])].append(tuple(add))

				self.__factorizations[str(pts[n])] = [list(x) for x in Set(self.__factorizations[str(pts[n])])]
				self.__max[str(pts[n])] = max(map(sum,self.__factorizations[str(pts[n])]))
			else:
				self.__factorizations[str(pts[n])] = []

	def Factorizations(self,v):
		if v in self.semigroup:
			if str(v) in self.__factorizations:
				return self.__factorizations[str(v)]
			else:
				genstr = ','.join(str(x) for x in self.gens)
				ret =  ConvertGapToSage(gap('FactorizationsVectorWRTList('+ str(v) +',['+ genstr + ']);'))
				self.__factorizations[str(v)] = ret
				return ret
		else:
			return []

	def LengthSet(self,v):
		return sorted(list(set(map(sum,self.Factorizations(v)))))
		
	def MinimalPresentation(self):
		return self.GenericGapCallGlobal('MinimalPresentationOfAffineSemigroup')

	def ElementsOfMinimalPresentation(self):
		m = self.MinimalPresentation()
		ret = []

		for p in m:
			x = y = 0
			for i in range(len(self.gens)):
				x += p[0][i]*self.gens[i][0]
				y += p[0][i]*self.gens[i][1]
			ret.append([x,y])

		return ret

	def OmegaPrimality(self, v=None):
		if v != None:
			return ConvertGapToSage(gap('OmegaPrimalityOfElementInAffineSemigroup('+ str(v) + ',' + self.semigroup.name() + ');'))
		else:
			return self.GenericGapCallGlobal('OmegaPrimalityOfAffineSemigroup')

	def DeltaSet(self,v):
		l = self.LengthSet(v);
		return [l[i+1]-l[i] for i in range(len(l)-1)]


	def Elasticity(self,v=None):
		if v!= None:
			len = self.LengthSet(v);
			return (max(len))/float(min(len))
		else:
			return self.GenericGapCallGlobal('ElasticityOfAffineSemigroup')

	def CatenaryDegree(self,v=None):
		if v!= None:
			return ConvertGapToSage(gap('CatenaryDegreeOfSetOfFactorizations('+ str(self.Factorizations(v)) +');'))
		else:
			return self.GenericGapCallGlobal('CatenaryDegreeOfAffineSemigroup')

	def __repr__(self):
		return '<Affine semigroup generated by %s>' % (self.gens)

	def BelongsToAffineSemigroup(self,v):
		if ConvertGapToSage(gap('BelongsToAffineSemigroup('+ str(v) + ',' + self.semigroup.name() + ');')) == True:
			return True
		else:
			return False