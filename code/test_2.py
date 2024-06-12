import numpy as np
from kp_to_point import kp_to_point

test_refarr = np.array([[323.311193, 354.019465, 0.0],
                        [528.986042, 653.911151, 362.668393],
                        [1011.141836, 835.867905, 876.634508],
                        [1139.267152, 1172.824856, 1236.164168],
                        [1476.439035, 1236.846677, 1578.443133],
                        [1810.239200, 960.541977, 2010.607450],
                        [2197.986866, 957.172407, 2397.336663],
                        [2184.499991, 589.889330, 2763.888932],
                        [2740.833599, 259.671518, 3409.124014],
                        [2875.702352, 57.497347, 3651.509728],
                        [2892.560946, -141.307254, 3850.498439],
                        [2858.843758, -340.111855, 4051.606909],
                        [2673.399222, -505.220761, 4299.243152],
                        [2391.965305, -1119.078801, 4972.745895],
                        [2391.965305, -1402.069982, 5254.984465],
                        [2503.916901, -1750.873066, 5620.339216],
                        [2247.086769, -2079.932578, 6036.651903],
                        [1654.401848, -2191.812813, 6638.196785],
                        [1002.448435, -2356.342569, 7308.793487],
                        [1035.375375, -2757.795175, 7710.516118],
                        [607.325155, -3027.623975, 8215.159248],
                        [139.762606, -3001.299214, 8682.205218],
                        [80.494114, -2626.171370, 9060.965601],
                        [-130.238302, -2382.667330, 9382.128122],
                        [-413.409987, -2073.351388, 9800.358799],
                        [-499.020031, -1750.873066, 10133.108314],
                        [-367.312271, -1507.369026, 10409.203819],
                        [-380.483047, -1040.104518, 10875.394832],
                        [-479.263867, -829.506430, 11107.381977]])

test_pointarr = np.array([[-89.711643, 256.474562],
                          [86.972542, 556.646209],
                          [469.788278, 556.646209],
                          [675.919827, 533.10333],
                          [858.493486, 562.531927],
                          [905.609269, 709.674891],
                          [876.161904, 903.903604],
                          [1035.177671, 1092.246597],
                          [1264.867113, 1074.589442],
                          [1364.988151, 1180.532376],
                          [1482.777608, 1292.361028],
                          [1694.798631, 1292.361028],
                          [1830.256507, 1045.160849],
                          [1895.040708, 774.417795],
                          [2018.719638, 709.674891],
                          [2324.972226, 715.560610],
                          [2560.551140, 621.389113],
                          [2743.124799, 374.188933],
                          [2701.898489, 109.331598],
                          [2872.693201, -214.382922],
                          [2949.256349, -420.383072],
                          [2937.477403, -761.754748],
                          [2778.461636, -861.811964],
                          [2360.309063, -961.869179],
                          [2236.630133, -1191.412203],
                          [2260.188025, -1403.298071],
                          [2631.224814, -1597.526783],
                          [2519.201075, -2011.279496],
                          [2401.411618, -2046.593808],
                          [2136.385340, -1993.622341],
                          [1730.011713, -2040.708089],
                          [1665.227511, -2181.965335],
                          [1612.222256, -2329.108299],
                          [1494.432799, -2452.708388],
                          [1093.948645, -2334.994017],
                          [905.485513, -2282.022550],
                          [840.701312, -2452.708388],
                          [828.922366, -2758.765753],
                          [876.038149, -3064.823118],
                          [569.785561, -3306.137579],
                          [328.317174, -3312.023298],
                          [186.969825, -3094.251711],
                          [228.196135, -2758.765753],
                          [27.954058, -2523.337011],
                          [-390.198514, -2358.536891],
                          [-472.651134, -2181.965335],
                          [-378.419568, -1981.850904],
                          [-284.188003, -1746.422161],
                          [-325.414313, -1616.936353],
                          [-560.993227, -1416.821922],
                          [-637.556374, -1293.221833],
                          [-431.424824, -1081.335964],
                          [-284.188003, -910.650126],
                          [-313.635367, -634.021354],
                          [-531.545863, -475.106953],
                          [-967.366854, -339.735426],
                          [-1450.303628, 166.436370],
                          [-1644.656232, 225.293555]])

maxoff_list = [100, 200, 500, 1000]


def test_without_maxoff(refarr, pointarr):
    print(f'Maxoff not specified: {kp_to_point.go(refarr, pointarr)}')


def test_with_maxoff(refarr, pointarr, maxoff):
    print(f'Maxoff = {maxoff}: {kp_to_point.go(refarr, pointarr, maxoff)}')

if __name__ == '__main__':
    for i, point in enumerate(test_pointarr):
        print(f'Point {i}')
        test_without_maxoff(test_refarr, np.array([point]))
        for maxoff in maxoff_list:
            test_with_maxoff(test_refarr, np.array([point]), maxoff)


    print('OK')


