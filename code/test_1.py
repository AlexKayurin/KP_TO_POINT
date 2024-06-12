import numpy as np
from kp_to_point import kp_to_point

test_refarr_1 = np.array([[1, -1, 0], [0, 0, 1.4142], [-1, 1, 2.8284]])

test_refarr_1_scaled = np.array([[1, -1, 0], [0, 0, 1], [-1, 1, 2]])

test_refarr_2 = np.array([[-1, -1, 0], [0, 0, 1.4142], [1, 1, 2.8284]])

test_refarr_1_scaled = np.array([[-1, -1, 0], [0, 0, 2], [1, 1, 4]])

test_pointarr = np.array([[-1, 1], [1, 1], [1, -1], [-1, -1], [0, 0], [-2, 2], [2, 2], [2, -2], [-2, -2]])

def test(refarr, pointarr):
    print(kp_to_point.go(refarr, pointarr))


if __name__ == '__main__':
    for refarr in [test_refarr_1, test_refarr_1_scaled, test_refarr_2, test_refarr_1_scaled]:
        test(refarr, test_pointarr)

    print('OK')


