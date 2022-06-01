from org.ejml.simple import SimpleMatrix

transformation_to_raw_registered = SimpleMatrix(3,4,True,[1.0, 0.0, 0.0, -7.0, 0.0, 1.0, 0.0, -7.0, 0.0, 0.0, 1.0, -46.0]);

rigid_transf = transformation_to_raw_registered.extractMatrix(0,3,0,3)
affine_transf = transformation_to_raw_registered.extractMatrix(0,3,3,4)

inverse_rigid = rigid_transf.invert()
inverse_affine = inverse_rigid.negative().mult(affine_transf)

full_inverse_transf = inverse_rigid.combine(0,inverse_rigid.numCols(),inverse_affine)


leg_point = SimpleMatrix(4,1,True, [160, 144, 0, 1])

new_coords = transformation_to_raw_registered.mult(leg_point)

print(transformation_to_raw_registered)
print(full_inverse_transf)

print(new_coords)

transf_to_uncropped_py = [[1.0, 0.0, 0.0, -7.0], [0.0, 1.0, 0.0, -7.0], [0.0, 0.0, 1.0, -46.0]]
leg_point_py = [[160], [144], [0]]

def point_coord_to_matrix(X, Y, Z):
	pass

def multiply_matrices(X, Y):
	C = [[0 for col in range(len(Y[0]))] for row in range(len(X))]
	print(C)
	# iterate through rows of X
	for i in range(len(X)):
		# iterate through columns of Y
		for j in range(len(Y[0])):
			# iterate through rows of Y
			for k in range(len(Y)):
				print(i, j)
				C[i][j] += X[i][k] * Y[k][j]

	return C

print(multiply_matrices(transf_to_uncropped_py, leg_point_py))



# A = A.transpose()
# for i in range(A.getNumElements()):
# 	print A.get(i)