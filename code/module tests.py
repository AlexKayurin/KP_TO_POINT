from kp_to_point import kp_to_point
from inspect import getmembers, isfunction, isclass


if __name__ == '__main__':
    print(getmembers(kp_to_point, isclass))
    print(getmembers(kp_to_point, isfunction))
