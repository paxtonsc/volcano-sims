import numpy as np

def compute_ref_mapping(x_tri):
	''' Compute jacobian
	x_tri: triangle nodes (..., 3, 2)

	Returns
	Matrix collection (..., 2, 2) such that reference coordinates can be
	obtained from

		M @ (x - x_tri[0,:]).

	'''
	_a00 = x_tri[...,1:2,0:1] - x_tri[...,0:1,0:1]
	_a01 = x_tri[...,2:3,0:1] - x_tri[...,0:1,0:1]
	_a10 = x_tri[...,1:2,1:2] - x_tri[...,0:1,1:2]
	_a11 = x_tri[...,2:3,1:2] - x_tri[...,0:1,1:2]
	_M = np.zeros((*x_tri.shape[:-2], 2, 2))
	_M[...,0:1,0:1] = _a11
	_M[...,0:1,1:2] = -_a01
	_M[...,1:2,0:1] = -_a10
	_M[...,1:2,1:2] = _a00
	_M /= (_a00 * _a11 - _a01 * _a10)
	return _M