

###### Math functions

def midpoint(x1, y1, x2, y2):
	x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)
	return ((x1 + x2) / 2, (y1 + y2) / 2)

def get_4_4_identity_matrix():
	identity_matrix = [[0 for col in range(4)] for row in range(4)]
	identity_matrix[0][0] = 1
	identity_matrix[1][1] = 1
	identity_matrix[2][2] = 1
	identity_matrix[3][3] = 1
	return identity_matrix

def multiply_matrices(X, Y):
	C = [[0 for col in range(len(Y[0]))] for row in range(len(X))]
	# iterate through rows of X
	for i in range(len(X)):
		# iterate through columns of Y
		for j in range(len(Y[0])):
			# iterate through rows of Y
			for k in range(len(Y)):
				C[i][j] += X[i][k] * Y[k][j]

	return C

def inverse_matrix(m):
	def return_transpose(mat):
		return map(list,zip(*mat))

	def return_matrix_minor(mat,i,j):
		return [row[:j] + row[j+1:] for row in (mat[:i]+mat[i+1:])]

	def return_determinant(m):
		if len(m) == 2:
			return m[0][0]*m[1][1]-m[0][1]*m[1][0]

		determinant = 0
		for c in range(len(m)):
			determinant += ((-1)**c)*m[0][c]*return_determinant(return_matrix_minor(m,0,c))
		return determinant
	
	determinant = return_determinant(m)
	if len(m) == 2:
		return [[m[1][1]/determinant, -1*m[0][1]/determinant],
				[-1*m[1][0]/determinant, m[0][0]/determinant]]

	cfs = []
	for r in range(len(m)):
		cfRow = []
		for c in range(len(m)):
			minor = return_matrix_minor(m,r,c)
			cfRow.append(((-1)**(r+c)) * return_determinant(minor))
		cfs.append(cfRow)
	cfs = return_transpose(cfs)
	for r in range(len(cfs)):
		for c in range(len(cfs)):
			cfs[r][c] = cfs[r][c]/determinant
	return cfs

