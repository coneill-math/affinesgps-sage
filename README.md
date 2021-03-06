# affinesgps-sage
Written by Madeline Chen and Christopher O'Neill.  

Provides Python classes for obtaining [cone decompositions](http://arxiv.org/abs/1503.08351) of quasipolynomial functions in the computer algebra system [Sage](http://sagemath.org/).  The package was developed as part of a senior thesis project by Madeline Chen, an undergraduate student at University of California Davis.  For a more detailed overview, see [Madeline's senior thesis](https://www.math.ucdavis.edu/files/3715/2903/7550/s18-chen-madeline-thesis.pdf).  

Please note that this is an *alpha version* and subject to change without notice.  

## License
affinesgps-sage is released under the terms of the [MIT license](https://tldrlegal.com/license/mit-license).  The MIT License is simple and easy to understand and it places almost no restrictions on what you can do with this software.

## Usage
The `AffineSemigroup` class supports the basic functionality of affine semigroups, and runs several core functionalities through GAP.

	sage: A = AffineSemigroup([[1,3],[3,3],[3,1]])
	sage: A.gens
	[[1, 3], [3, 1], [3, 3]]
	sage: [6,6] in A
	True
	sage: [4,3] in A
	False
	sage: A.Plot(20,20)
	Launched png viewer for Graphics object consisting of 1 graphics primitive
	sage: A.Plot(20,20,MaxFact)
	Launched png viewer for Graphics object consisting of 182 graphics primitives
	sage: A.Factorizations([40,50])
	[[13, 8, 1], [10, 5, 5], [7, 2, 9]]

The `FactorizationsUpToElement()` member function utilizes a dynamic programming algorithm to generate the factorizations of every element up to the element.  For each element, this function saves the factorization of the each element interally.

	sage: A = AffineSemigroup([[1,3],[3,3],[3,1]])
	sage: time A.Factorizations([50,50])
	CPU times: user 33.1 ms, sys: 6.31 ms, total: 39.4 ms
	Wall time: 41.5 ms
	[[11, 11, 2], [8, 8, 6], [5, 5, 10], [2, 2, 14]]
	sage: time A.FactorizationsUpToElement([50,50])
	CPU times: user 3.58 s, sys: 599 ms, total: 4.18 s
	Wall time: 4.66 s
	sage: time A.Factorizations([50,50])
	CPU times: user 1.2 ms, sys: 703 µs, total: 1.9 ms
	Wall time: 1.46 ms
	[[11, 11, 2], [8, 8, 6], [5, 5, 10], [2, 2, 14]]

The `ConeDecomposition` class provides an automated process for decomposing quasipolynomial functions defined on affine semigroups.  Using the member function `Decomposition()`, the user can recursively define cones in terms of basepoints and generators. 

	sage: disjoint =  ConeDecomposition(A,30,30)
	sage: conelist = disjoint.Decomposition(30,30)
	Beginning Disjoint Decomposition for Affine Semigroup generated by [[1, 3], [3, 1], [3, 3]]
	Minimal Presentation = [[[3, 3, 0], [0, 0, 4]]] @ [[12, 12]]

At each iteration of defining a cone, the function gives a list of suggested generators, using the original generators of the affine semigroup. Cone definition is done by running the cone object previously defined. Furthermore, the function gives suggested basepoints, based off of the remaining points in the function, utilizing the \code{RemainingPts} function. 

	Defining Cone 1:
		Cone generators Options-  0:[1, 3]  1:[3, 1]  2:[3, 3]   
		Please enter generators:[0,1]
		Set Basepoint: (enter vector or 1 to keep set @ [0, 0])  1
		Preview -- Cone generate by [[1, 3], [3, 1]] with basepoint @ [0, 0]
	Launched png viewer for Graphics object consisting of 180 graphics primitives
		*** Warning -- Cone Overlap @ []
		0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program
		Choose option: 1
	
	Defining Cone 1:
		Cone generators Options-  0:[1, 3]  1:[3, 1]  2:[3, 3]   
		Please enter generators:[0,1]
		Set Basepoint: (enter vector or 1 to keep set @ [0, 0])  [3,3]
		Preview -- Cone generate by [[1, 3], [3, 1]] with basepoint @ [3, 3]
	Launched png viewer for Graphics object consisting of 180 graphics primitives
		*** Warning -- Cone Overlap @ []
		0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program
		Choose option: 3
		Finished with 1 cone.


The `Decomposition()` function also supports the input of unique functions.

	sage: maxFact = ConeDecomposition(A,15,15,MaxFact)
	sage: conelist =maxFact.Decomposition(15,15)
	Begining Decomposition for Affine Semigroup generate by [[1, 3], [3, 1], [3, 3]]
	Minimal Presentation = [[[3, 3, 0], [0, 0, 4]]] @ [[12, 12]]
	
	Defining Cone 1:
		Cone generators Options-  0:[1, 3]  1:[3, 1]  2:[3, 3]   
		Please enter generators:[0,2]
		Set Basepoint: (enter vector or 1 to keep set @ [0, 0])  1
		Preview -- Cone generate by [[1, 3], [3, 3]] with basepoint @ [0, 0]
	Launched png viewer for Graphics object consisting of 204 graphics primitives
		0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program
		Choose option: 0
		Cone generate by [[1, 3], [3, 3]] with basepoint @ [0, 0]

With each iteration of defining a new cone, the user has the option to: 0) keep the current cone, 1) remove the current cone, 2) redefine the current cone, 3) add the current cone to the list and end the program, 4) end the program. By providing a plot preview of the current decomposition process, this gives the user a visual insight to where they should place the next cone. In addition, with the given options at each cone iteration, this interface gives the user flexibility to redefine previous cones. 

	Defining Cone 2:
		Cone generators Options-  0:[1, 3]  1:[3, 1]  2:[3, 3]   
		Please enter generators:[0,1]
		Set Basepoint: (enter vector or 1 to keep set @ [3, 1])  1
		Preview -- Cone generate by [[1, 3], [3, 1]] with basepoint @ [3, 1]
	Launched png viewer for Graphics object consisting of 148 graphics primitives
		*** Warning -- Cone Overlap @ [[12, 12], [13, 15], [14, 18], [15, 21], [16, 24], [17, 27], [18, 30], [24, 24], [25, 27], [26, 30]]
		0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program
		Choose option: 1
	
	Defining Cone 2:
		Cone generators Options-  0:[1, 3]  1:[3, 1]  2:[3, 3]   
		Please enter generators:[1,2]
		Set Basepoint: (enter vector or 1 to keep set @ [3, 1])  1
		Preview -- Cone generate by [[3, 1], [3, 3]] with basepoint @ [3, 1]
	Launched png viewer for Graphics object consisting of 150 graphics primitives
		0:Keep Cone / 1:Remove Current Cone / 2:Redefine Previous Cone(s) / 3: Add Cone & End Program / 4: End Program
		Choose option: 0

You can find action shots in the `images` folder.  
