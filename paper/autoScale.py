## \file   autoScale.py
#  \brief  a prograom for automatic finite size scaling analyses
#
#	autoScale.py -- a program for automatic finite size scaling analyses
#	Copyright (C) 2007-2009  Oliver Melchert
#
#	This program is free software; you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation; either version 2 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program; if not, write to the Free Software
#	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
#
#	On Debian systems, the complete text of the GNU General
#	Public License can be found in `/usr/share/common-licenses/GPL'
#
#  \date   02.10.09
#  \author Oliver Melchert
import sys,os,math

pow=math.pow
sqrt=math.sqrt
def div(a,b): return(float(a)/float(b))


def rootBisection(myFunc,xMin,xMax,epsilon):
	"""shrink interval that contains root of function, method
	does not fail but is quite slow

	NOTE: used for root search in errorAnalysis()

	Input:
	myFunc			-- objective function 
	xMin,xMax		-- interval that contains root
	epsilon			-- extension of minimal interval that
				   signals convergence

	Returns: (xMin,xMax)
	xMin,xMax		-- interval with small extension that
				   contaings root
	"""
	# initialize function value at interval boundaries
	fMin = myFunc(xMin)
	fMax = myFunc(xMax)
	#file=sys.stdout
	#file.write("# epsilon=%lf\n"%(epsilon))
	#file.write("# iter: [xMin,xMax]; [fMin,fMax]\n")

	iter=0
	# as long as distance of best function value
	# encountered so far is further away from the
	# desired function value than epsilon: iterate 
	while(abs(xMax-xMin)>2.*epsilon):
	  midPt  = div(xMin+xMax,2)	
	  fMidPt = myFunc(midPt)	
	   
	  # if "left" half interval does not contain
	  # the root: discard that part of the interval
	  if(fMin*fMidPt > 0.): xMin=midPt; fMin=fMidPt	
	  # else: root is contained in "left" half of
	  # interval, hence discard right half 
	  else: xMax=midPt; fMax=fMidPt
	 
	  iter+=1
	  #file.write("%3d: [%10.9lf,%10.9lf] [%10.9lf,%10.9lf]\n"%(iter,xMin,xMax,fMin,fMax))

	# return final value xMin +- epsilon 
	return div(xMin+xMax,2)

def getBrackets(dumFunc,midVal):
	"""find bracketing intervals of roots 

	NOTE: there are two roots and val lies in between them.
	So, this is a simple routine to find the intervals that
	bracket both roots

	Input:
	dumFunc		-- function that changes sign at roots
	val		-- value that is located in betwee roots

	Returns: (lBrack,rBrack)
	lBrack,rBrack	-- left/right interval containing root
	"""
	fac=0.01 # factor used to modify value of midVal
			 # so as to bracket an interval that 
			 # contains a root of the objective function
	
	midF = dumFunc(midVal)
	lVal=midVal*(1.-fac)
	lF  =dumFunc(lVal)
	while(lF*midF>0.):
		lVal *= (1.-fac)
		midF=lF
		lF=dumFunc(lVal)
		print "lVal, lF", lVal , lF 
	lBrack=[lVal,midVal]

	midF = dumFunc(midVal)
	rVal=midVal *(1.+fac)
	rF  =dumFunc(rVal)
	while(rF*midF>0.):
		rVal *= (1.+fac)
		midF=rF
		rF=dumFunc(rVal)
		print "rVal, rF", rVal, rF 
	rBrack=[midVal,rVal]

	return lBrack,rBrack

def amoeba(my_func,p,y,opt_flag,ftol,report_fitness):
	''' Simplex algorithm of Nelder & Mead, adopted from Numerical recipes in C.
	
	NOTE: multidimensional minimization of function my_func(q), where q
	is ndim (number of scaling parameters) vector.

	Input:
	p[ndim+1][ndim] 	-- matrix of ndim+1 initial vertices of simplex
	y[ndim+1] 	  	-- fitness of the ndim+1 starting vertices
	opt_flag[ndim]  	-- list of optimization parameters
	  		    	   opt_flag[par] = 1 (0) -> (don't) optimize par
	ftol		  	-- stopping tolerance of simplex extension
	report_fitness 		-- write best fitness at each step to stdout	

	Returns: (p,y)
	p[][] and y[] 		-- new vertices within ftol of minimal fitness value
	  		   	   the best result is reported in 0th-slot p[0][],y[0]
	'''
	# set parameters
	TINY = 1.0e-10	# simply a small number
	ITMAX = 800	# maximum allowed iterations
	
	ndim = len(y)-1 # number of scaling parameters
	iter=0		# start iterations counter
	
	psum = [1.0 for i in range(ndim)] # create psum array
	
	# get psum
	for n in range(ndim):
		sum = 0.
		for m in range(ndim+1):
			sum += p[m][n]
		psum[n]=sum
		
	# start iteration
	while(1):
	 	# get worst (ihi), next worst (inhi) 
		# and best (ilo) vertex	in simplex
		ilo=1
		if(y[0] > y[1]):
			ihi =0
			inhi=1
		else:
			ihi =1
			inhi=0
		
		for i in range(ndim+1):
			if(y[i] < y[ilo]) : ilo = i
			if(y[i] > y[ihi]) : 
				inhi = ihi
				ihi  = i
			elif( (y[i] > y[inhi]) and (i != ihi) ): inhi = i

		# determine fractional range from highest to lowest value
		rtol = 2.0 * abs(y[ihi]-y[ilo]) / 1.0*(abs(y[ihi]) + abs(y[ilo]) + TINY)
		# report best fitness in current simplex to stdout
		if(report_fitness): print iter, y[ilo], rtol	
		# return if satisfactory, 
		# first put best vertex to slot 0
		if(rtol < ftol): 
			y[0],y[ilo] = y[ilo],y[0]
			for n in range(ndim):
				p[0][n],p[ilo][n]   = p[ilo][n],p[0][n]
			
			return p[0],y[0],iter
		
		# report if no iterations left 
		if(iter > ITMAX):  
			print "ITMAX exceeded in amoeba" 
			sys.exit(1)
		
		iter += 2	# start new iteration
		
		# reflect simplex from the high point
		# extrapolate with factor -1 through face
		ytry = amotry(my_func,p,y,psum,ndim,ihi,-1.0,opt_flag)
		if(ytry <= y[ilo]):	
			# result superceeds best vertex
			# try another extrapolation by factor 2
			ytry = amotry(my_func,p,y,psum,ndim,ihi,2.0,opt_flag)
		elif(ytry >= y[inhi]):
			# reflected point is worse than second highest
			# look for intermediate low point (1-d contraction)
			ysave = y[ihi]
			ytry = amotry(my_func,p,y,psum,ndim,ihi,0.5,opt_flag)
			if(ytry >= ysave):
				# no improvement ... 
				# rather contract arround best point
				for i in range(ndim+1):
					if(i != ilo):
						for j in range(ndim):
						  # don't change the p[i][j] where opt_flag[j] == 0
						  if opt_flag[j]!=0: # changed mo, 22.01
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j])
						  else: # mel
							psum[j]=p[i][j] #mel

						y[i]=my_func(psum) 
				iter+=1			
				# get psum
				for n in range(ndim):
					sum = 0.
					for m in range(ndim+1):
						sum += p[m][n]
					psum[n]=sum
		else: 
			iter -=1

def amotry(my_func,p,y,psum,ndim,ihi,fac,opt_flag):
	'''Extrapolate by factor fac through simplex face
	takes worst point and tries to improve. upon success
	the bad point is replaced by the good point
	
	Input:
	my_func			-- quality function that measures fitness
	p			-- list of points that define simplex
	y			-- fitness values of simplex points
	psum			-- particular point of the simplex
	ndim			-- number of scaling parameters
	ihi			-- index of best simplex point
	fac			-- scaling factor for simplex manipulation
	opt_flag		-- tells if scaling parameter will be optimized

	Returns: (ytry)
	ytry			-- fitness of manipulated simplex point
	'''
	fac1 = (1.0-fac)/(1.0*ndim)
	fac2 = fac1 -fac
	ptry = [1.0 for i in range(ndim)] # create ptry array
	for j in range(ndim):
	   # don't change the ptry[j], where opt_flag[j]==0
	   if opt_flag[j]!=0:
		ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2
	   else:
		ptry[j] = p[ihi][j]
	ytry = my_func(ptry)
	if(ytry<y[ihi]):
		y[ihi]=ytry
		for j in range(ndim):
			psum[j] += ptry[j]-p[ihi][j]
			p[ihi][j] = ptry[j]
	
	return ytry

def iniSimplex(myFunc,trialPoint,delta,optFlag):
	'''Construct initial simplex

	Input:
	myFunc		-- function that implements scaling assumption
	trialPoint	-- initial guess as trial point for simplex 
			   construction
	delta		-- scaling parameter for simplex points
	optFlag		-- indicates if a scaling parameter will be
			   optimized

	Returns: (p,y)
	p		-- initial simplex
	y		-- function values of simplex points
	'''
	# create dim+1 x dim array
	dim = len(trialPoint)
	p = [ [1.0 for i in range(dim)] for j in range(dim+1) ]
	y = [1.0 for i in range(dim+1)]
	for i in range(dim+1):
		# construct vertices that define an initial simplex
		dum=0.
		for j in range(dim):
			if( (i-1) == j ) and (optFlag[j]!=0): 
				dum = trialPoint[j]*(1+delta)
			else : 	dum = trialPoint[j]
			# fill vertex matrix
			p[i][j] = dum
		# get function values belonging to vertices
		y[i] = myFunc(p[i])		
	return p,y


class myValue:
	"""Returns new instance of class 'myValue' and assigns 
	the respective object to the local variable"""
	def __init__(self,L,x,y,dy):
 		"""Define initial state for a new created instance 
		of the 'myValue' class"""
		self.L  = L
		self.x  = x
		self.y  = y
		self.dy = dy

	def __repr__(self):
		return '%lf %lf %lf %lf'%(self.L,self.x,self.y,self.dy)

class myRawData:
	"""Returns new instance of class 'myData' and assigns 
	the respective object to the local variable
	
	Uses:
	class myValue		-- used in method fetchData()
	"""
	def __init__(self,data):
 		"""Define initial state for a new created instance 
		of the 'myData' class"""
		self.dataSet = data

	def listDataSets(self):
		"""list data sets"""
		for L,data in self.dataSet.iteritems():
		  for val in data:
		    print L,val.x,val.y,val.dy

	def listDataSetsScaled(self,scaleAssumption):
		"""list scaled data sets"""
		for L,data in self.dataSet.iteritems():
		  for val in data:
		    print scaleAssumption.scale(val)


class myScaleAssumption:
	"""Returns new instance of class 'myScaleAssumption' and assigns 
	the respective object to the local variable

	Uses:
	myValue			-- data structure for data point
	"""
	def __init__(self, xc, a, b, c, w, scaledXMin, scaledXMax, xco=True, ao=True, bo=True, co=True, wo=True):
 		"""Define initial state for a new created instance 
		of the 'myScaleAssumption' class"""
		# scaling parameters/corresp. optimization flags
		self.xc = xc; self.xco = xco
		self.a  = a; self.ao  = ao
		self.b  = b; self.bo  = bo
		self.c  = c; self.co  = co
		self.w  = w; self.wo  = wo
		# intervall on rescaled x-axes	
		self.scaledXMin = scaledXMin
		self.scaledXMax = scaledXMax

	def updateByHand(self,xc,a,b):
		"""update scaling parameters
		
		Input:
		xc		-- critical point
		a,b		-- crit exponents

		Returns: nothing
		"""
		self.xc = xc
		self.a  = a
		self.b  = b
		self.c  = c
		self.w  = w
		
	def updateFromList(self,par):
		"""update scaling parameters from supplied
		parameter list
		
		Input:
		par		-- list of scaling parameters, must
				   take form par=[xc,a,b]

		Returns: nothing
		"""
		self.xc = par[0]
		self.a  = par[1]
		self.b  = par[2]
		self.c  = par[3]
		self.w  = par[4]
		
	def scaleParNames(self):
		"""return list of scaling parameter Names

		Input: nothing

		Returns: (nameList)
		nameList	-- list of scaling parameters
		"""
		return ['xc','a','b','c','w']
		
	def scaleParList(self):
		"""return list of scaling parameters and optimization flags

		Input: nothing

		Returns: (parList,optFlag)
		parList		-- list of scaling parameters
		optFlag		-- corresponding optimization flag
		"""
		return [self.xc,self.a,self.b, self.c, self.w],[self.xco,self.ao,self.bo, self.co, self.wo]

		
	def scale(self,val):
		"""apply scaling assumption 

		x  => (x-xc) L^a
		y  =>  y L^b / (1+c*L^(-w))
		dy => dy L^b / (1+c*L^(-w))

		to supplied data point (L,x,y,dy)
		
		Input:
		val		-- data structure representing unscaled data point

		Returns: (scaledVal)
		scaledVal	-- input value after scaling assumption was applied
		"""
		return myValue(val.L, \
                (val.x-self.xc)*pow(val.L,self.a),\
			    val.y*pow(val.L,self.b)/(1.+self.c*val.L**(-self.w)),\
			    val.dy*pow(val.L,self.b)/(1.+self.c*val.L**(-self.w)))

	def listScalePar(self,fileStream=sys.stdout):
		"""list scaling parameters

		Input:
		fileStream		-- file out stream

		Returns: nothing, but writes scaling parameters to
		         specified file out stream
		"""
		fileStream.write("xc=%f a=%f b=%f c=%f w=%f\n"%(self.xc,self.a,self.b, self.c, self.w))
	

class myFunc(myRawData,myScaleAssumption):
	"""Returns new instance of class 'myFunc' and assigns 
	the respective object to the local variable
	
	Inherits:
	class myRawData		-- implements raw data container 
	class myScaleAssumption	-- implements scaling assumption
	"""
	def __init__(self, data, xc, a, b, c, w, scaledXMin, scaledXMax, xco=True, ao=True, bo=True, co=True, wo=True):
 		"""Define initial state for a new created instance 
		of the 'myFunc' class. Takes care that base classes
		are initialized properly
		"""
		myRawData.__init__(self, data)		# ini instance of class myRawData
		myScaleAssumption.__init__(self, xc, a, b, c, w, scaledXMin, scaledXMax, xco, ao, bo, co, wo)	# ini instance of class myScale Assumption

		
	def scaleData(self,scalePar):
		"""scale data sets so as to compute quality function

		Input:
		scalePar		-- list of scaling parameters in the form 
					   required by method updateFromList() of
					   base class myScaleAssumption

		Returns: (S)
		S			-- quality of the data collapse given
					   the scaling assumption in myValue.scaled()
		"""
		SList = []
		self.updateFromList(scalePar)
		for L,data in self.dataSet.iteritems():
		  for val in data:
		     # scaled value
		     sVal = self.scale(val)
		     x,y,dy = sVal.x,sVal.y,sVal.dy
		     # restrict scaling analysis on rescaled abscissa
		     if self.scaledXMin <= x <= self.scaledXMax:
			# select subset for estimation of master curve at x
			subSet = self.selectSubset(val)	
			if subSet != []:
			  	# linear regression to estimate master curve at x
				Y,dY2 = self.llsFit(x,subSet)
				# extimate 'quality' of the current point
				chi2 = div((y-Y)*(y-Y),dy*dy+dY2)
				SList.append(float(chi2))
			del sVal
			del subSet
			
		# determine quality of the data collapse
		if len(SList)>0: 
			S=div(sum(SList),len(SList))
		else:		 
			S=99999.99

		del SList
		return S
	
	def selectSubset(self,pivVal):
		"""select subset of data points (x,y,dy) that serve to 
		perform linear regression in order to compute value Y,dY of 
		the unknown master curve at the pivoting point pivVal.x

		NOTE: x-values contained in raw data sets do not need to
		be supplied in increasing order

		Input:
		pivVal			-- pivoting value for which the value
					   of the unknown master curve is to 
					   be computed

		Returns: (subSet)
		subSet			-- subset of data triples (x,y,dy) for
					   linear regression 
		"""
		# scale pivoting value
		sPiv = self.scale(pivVal)
		# initialize empty subset 
		subSet=[]
		for L in [L for L in self.dataSet.keys() if L!=pivVal.L]:
		  # initialize 'empty' min/max pair for this system size
		  maxVal,minVal=None,None
		  for val in self.dataSet[L]:
			sVal=self.scale(val)
			# get largest smaller value
			if sVal.x<=sPiv.x and sVal.x>=minVal:
			   minVal=sVal
			# get smallest larger value
			if sVal.x>sPiv.x and (maxVal==None or sVal.x<=maxVal):
			   maxVal=sVal
			#elif sVal.x>sPiv.x and sVal.x<=maxVal:
			#   maxVal=sVal
			del sVal
		  # add min/max pair to current set if both 'exist'
		  if minVal!=None and maxVal!=None:
			  subSet+=[minVal,maxVal]
		del sPiv
		return subSet
	
		
	def llsFit(self,pivX,subSet):
		"""perform linear least squares fit of the data points
		in subSet to a straight line to result in estimate Y=A+B*x 

		Input:
		pivX			-- x value for which Y,dY will be computed
		subSet			-- list of data points (x,y,dy)

		Return: (Y,dY2)
		Y,dY2			-- estimate of master function with error
		"""
		K=Kx=Ky=Kxx=Kxy=0.
		for val in subSet:
			x,y,dy=val.x,val.y,val.dy
			wgt = div( 1.,dy*dy)
			K   += wgt 
			Kx  += x*wgt
			Ky  += y*wgt
			Kxx += x*x*wgt
			Kxy += x*y*wgt
		
		fac = K*Kxx -Kx*Kx
		A   = div( Ky*Kxx - Kx*Kxy  ,fac)
		B   = div( K*Kxy-Kx*Ky  ,fac)
		Y=A+B*pivX
		dY2 = div( Kxx-2.*pivX*Kx+pivX*pivX*K,fac    )
		return Y,dY2 

	def listState(self,scalePar,fileStream=sys.stdout):
		"""list scaling parameters

		Input:
		fileStream		-- file out stream

		Returns: nothing, but writes scaling parameters to
		         specified file out stream
		"""
		quality=self.scaleData(scalePar)
		fileStream.write("dx = [%lf:%lf]  xc = %f  a = %f  b = %f  c = %f w = %f S = %f\n"%\
			(self.scaledXMin,self.scaledXMax,self.xc,self.a,self.b,self,c, self.w, quality))

		
def errorAnalysis(f):
	"""perform error analysis for scaling parameters

	Input: 
	f		-- data structure containing raw data and scaling assumption

	Returns: (parErr)
	parErr		-- errors for scaling parameters
	"""
	# get list of scaling parameters and corresponding optimization flags
	scalePar,optFlag=f.scaleParList()
	bestS = f.scaleData(scalePar)
	# save error bounds for scaling parameters
	parErr = {}
	# print msg
	print "# PERFORM S+1 ERROR ANALYSIS FOR SCALING PARAMETERS"
	
	for parId in range(len(optFlag)):
	  # if scaling parameter was optimized determine error bounds
	  if optFlag[parId]==1:
		# define dummy function that serves to bracket minimum 
		# for S+1 error analysis
		def dumFunc(val): 
			scalePar[parId]=val
			return f.scaleData(scalePar)-(bestS+1.)

		# save original value of the scaling paramter
		pivVal  = scalePar[parId]	

		# find bracketing intervalls of roots
		lBrack,rBrack=getBrackets(dumFunc,pivVal) ; print 'bracket done', parId, len(optFlag)

		# estimate parameter value at root
		err1 = rootBisection(dumFunc,lBrack[0],lBrack[1],10**(-5))	; print 'error1 done', parId, len(optFlag)
		scalePar[parId]=pivVal		# restore original value

		# estimate parameter value at root
		err2 = rootBisection(dumFunc,rBrack[0],rBrack[1],10**(-5)) ; print 'error2 done', parId, len(optFlag)
		scalePar[parId]=pivVal		# restore original value

		parErr[parId]=[scalePar[parId],abs(pivVal-err1),abs(pivVal-err2)]

	  # if scaling parameter was not optimized return error 0.
	  else:
		parErr[parId]=[scalePar[parId],0.,0.]

	return parErr
	

def autoScale(data, xc, a, b, c, w, scaledXMin, scaledXMax, xco=True, ao=True, bo=True, co=True, wo=True, repFit=1, getError=1):
	
	# set default parameters for minimization routine 
	outFile	 = sys.stdout	# default outfile
	delta	 = 0.001		# scaling parameter for construction of initial simplex
	simpTol	 = 10**(-6)		# simplex tolerance: if simplex extension gets smaller
	        				# than this, minimization routine has 'converged'
	
	# initialize instance of myFunc that carries raw data 
	# and scaling assumption, must be initialized before 
	# command line parameters are parsed
	f = myFunc(data, xc, a, b, c, w, scaledXMin, scaledXMax, xco, ao, bo, co, wo)
  	
	# abbreviation for call to quality function
	S=f.scaleData
	########## END: PARSE COMMAND LINE ARGUMENTS #########################

	# get list of scaling parameters and corresponding optimization flags
	scalePar,optFlag=f.scaleParList()

	########## SCALING PARAMETER OPTIMIZATION #########################
	# initialize simplex
	p,y = iniSimplex(S,scalePar,delta,optFlag)
	# minimize quality function S by adjusting scaling parameters	
	pBest,yBest,nIter=amoeba(S,p,y,optFlag,simpTol,repFit)
	print pBest , yBest, nIter
	########## END: SCALING PARAMETER OPTIMIZATION #####################
	xc, a, b, c, w, S =  f.xc, f.a, f.b, f.c, f.w, S(pBest)

	if getError:
		# here, f must containt best scaling parameters
		parErr=errorAnalysis(f); print parErr  
		# list scaling parameters with associated +- errors
		scaleParNames=f.scaleParNames()
		outFile.write("# S+1 error analysis yields:\n")
		outFile.write("# Scaling analysis restricted to\n")
		outFile.write("%4s = [%lf : %lf]\n"%('xr',f.scaledXMin,f.scaledXMax))
		outFile.write("# <scalePar>  <-Err>  <+Err>\n")
		for parId,values in parErr.iteritems():
			outFile.write("%4s = %lf %lf %lf\n"%\
					(scaleParNames[parId],values[0],values[1],values[2]))
			
	else:
		f.listState(pBest,outFile)

	return xc, a, b, c, w, S 
