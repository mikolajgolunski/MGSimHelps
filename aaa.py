class Aaa:
    def __init__(self):
        self.bins = [
            [
                [
                    [1], [2], [3]
                ],
                [
                    [4], [5], [6]
                ],
                [
                    [7], [8], [9]
                ]
            ],
            [
                [
                    [10], [11], [12]
                ],
                [
                    [13], [14], [15]
                ],
                [
                    [16], [17], [18]
                ]
            ],
            [
                [
                    [19], [20], [21]
                ],
                [
                    [22], [23], [24]
                ],
                [
                    [25], [26], [27]
                ]
            ]
        ]

    def aaa(self):
        atomCounter = 14
        xmax = len(self.bins) - 1
        ymax = len(self.bins[0]) - 1
        zmax = len(self.bins[0][0]) - 1
        big_bin = []
        for x, binx in enumerate(self.bins):
            for y, biny in enumerate(binx):
                for z, binz in enumerate(biny):
                    for i, binatom in enumerate(binz):
                        if binatom == atomCounter:
                            big_bin.extend(self.bins[x][y][z])
                            if 0 < x < xmax:
                                big_bin.extend(self.bins[x-1][y][z])
                                big_bin.extend(self.bins[x+1][y][z])
                                if 0 < y < ymax:
                                    big_bin.extend(self.bins[x][y-1][z])
                                    big_bin.extend(self.bins[x][y+1][z])
                                    big_bin.extend(self.bins[x-1][y-1][z])
                                    big_bin.extend(self.bins[x+1][y-1][z])
                                    big_bin.extend(self.bins[x-1][y+1][z])
                                    big_bin.extend(self.bins[x+1][y+1][z])
                                    if 0 < z < zmax:
                                        big_bin.extend(self.bins[x][y][z-1])
                                        big_bin.extend(self.bins[x][y][z+1])
                                        big_bin.extend(self.bins[x-1][y][z-1])
                                        big_bin.extend(self.bins[x+1][y][z-1])
                                        big_bin.extend(self.bins[x-1][y][z+1])
                                        big_bin.extend(self.bins[x+1][y][z+1])
                                        big_bin.extend(self.bins[x][y-1][z-1])
                                        big_bin.extend(self.bins[x][y+1][z-1])
                                        big_bin.extend(self.bins[x-1][y-1][z-1])
                                        big_bin.extend(self.bins[x+1][y-1][z-1])
                                        big_bin.extend(self.bins[x-1][y+1][z-1])
                                        big_bin.extend(self.bins[x+1][y+1][z-1])
                                        big_bin.extend(self.bins[x][y-1][z+1])
                                        big_bin.extend(self.bins[x][y+1][z+1])
                                        big_bin.extend(self.bins[x-1][y-1][z+1])
                                        big_bin.extend(self.bins[x+1][y-1][z+1])
                                        big_bin.extend(self.bins[x-1][y+1][z+1])
                                        big_bin.extend(self.bins[x+1][y+1][z+1])
                                    else:
                                        if z > 0:
                                            big_bin.extend(self.bins[x][y][z-1])
                                            big_bin.extend(self.bins[x-1][y][z-1])
                                            big_bin.extend(self.bins[x+1][y][z-1])
                                            big_bin.extend(self.bins[x][y-1][z-1])
                                            big_bin.extend(self.bins[x][y+1][z-1])
                                            big_bin.extend(self.bins[x-1][y-1][z-1])
                                            big_bin.extend(self.bins[x+1][y-1][z-1])
                                            big_bin.extend(self.bins[x-1][y+1][z-1])
                                            big_bin.extend(self.bins[x+1][y+1][z-1])
                                        if z < zmax:
                                            big_bin.extend(self.bins[x][y][z+1])
                                            big_bin.extend(self.bins[x-1][y][z+1])
                                            big_bin.extend(self.bins[x+1][y][z+1])
                                            big_bin.extend(self.bins[x][y-1][z+1])
                                            big_bin.extend(self.bins[x][y+1][z+1])
                                            big_bin.extend(self.bins[x-1][y-1][z+1])
                                            big_bin.extend(self.bins[x+1][y-1][z+1])
                                            big_bin.extend(self.bins[x-1][y+1][z+1])
                                            big_bin.extend(self.bins[x+1][y+1][z+1])
                                else:
                                    if y > 0:
                                        big_bin.extend(self.bins[x][y-1][z])
                                        big_bin.extend(self.bins[x-1][y-1][z])
                                        big_bin.extend(self.bins[x+1][y-1][z])
                                        if 0 < z < zmax:
                                            big_bin.extend(self.bins[x][y][z-1])
                                            big_bin.extend(self.bins[x][y][z+1])
                                            big_bin.extend(self.bins[x-1][y][z-1])
                                            big_bin.extend(self.bins[x+1][y][z-1])
                                            big_bin.extend(self.bins[x-1][y][z+1])
                                            big_bin.extend(self.bins[x+1][y][z+1])
                                            big_bin.extend(self.bins[x][y-1][z-1])
                                            big_bin.extend(self.bins[x-1][y-1][z-1])
                                            big_bin.extend(self.bins[x+1][y-1][z-1])
                                            big_bin.extend(self.bins[x][y-1][z+1])
                                            big_bin.extend(self.bins[x-1][y-1][z+1])
                                            big_bin.extend(self.bins[x+1][y-1][z+1])
                                        else:
                                            if z > 0:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x-1][y][z-1])
                                                big_bin.extend(self.bins[x+1][y][z-1])
                                                big_bin.extend(self.bins[x][y-1][z-1])
                                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                            if z < zmax:
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x-1][y][z+1])
                                                big_bin.extend(self.bins[x+1][y][z+1])
                                                big_bin.extend(self.bins[x][y-1][z+1])
                                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                                big_bin.extend(self.bins[x+1][y-1][z+1])
                                    if y < ymax:
                                        big_bin.extend(self.bins[x][y+1][z])
                                        big_bin.extend(self.bins[x-1][y+1][z])
                                        big_bin.extend(self.bins[x+1][y+1][z])
                                        if 0 < z < zmax:
                                            big_bin.extend(self.bins[x][y][z-1])
                                            big_bin.extend(self.bins[x][y][z+1])
                                            big_bin.extend(self.bins[x-1][y][z-1])
                                            big_bin.extend(self.bins[x+1][y][z-1])
                                            big_bin.extend(self.bins[x-1][y][z+1])
                                            big_bin.extend(self.bins[x+1][y][z+1])
                                            big_bin.extend(self.bins[x][y+1][z-1])
                                            big_bin.extend(self.bins[x-1][y+1][z-1])
                                            big_bin.extend(self.bins[x+1][y+1][z-1])
                                            big_bin.extend(self.bins[x][y+1][z+1])
                                            big_bin.extend(self.bins[x-1][y+1][z+1])
                                            big_bin.extend(self.bins[x+1][y+1][z+1])
                                        else:
                                            if z > 0:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x-1][y][z-1])
                                                big_bin.extend(self.bins[x+1][y][z-1])
                                                big_bin.extend(self.bins[x][y+1][z-1])
                                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                                big_bin.extend(self.bins[x+1][y+1][z-1])
                                            if z < zmax:
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x-1][y][z+1])
                                                big_bin.extend(self.bins[x+1][y][z+1])
                                                big_bin.extend(self.bins[x][y+1][z+1])
                                                big_bin.extend(self.bins[x-1][y+1][z+1])
                                                big_bin.extend(self.bins[x+1][y+1][z+1])
                            else:
                                if x > 0:
                                    big_bin.extend(self.bins[x-1][y][z])
                                    if 0 < y < ymax:
                                        big_bin.extend(self.bins[x][y-1][z])
                                        big_bin.extend(self.bins[x][y+1][z])
                                        big_bin.extend(self.bins[x-1][y-1][z])
                                        big_bin.extend(self.bins[x-1][y+1][z])
                                        if 0 < z < zmax:
                                            big_bin.extend(self.bins[x][y][z-1])
                                            big_bin.extend(self.bins[x][y][z+1])
                                            big_bin.extend(self.bins[x-1][y][z-1])
                                            big_bin.extend(self.bins[x-1][y][z+1])
                                            big_bin.extend(self.bins[x][y-1][z-1])
                                            big_bin.extend(self.bins[x][y+1][z-1])
                                            big_bin.extend(self.bins[x-1][y-1][z-1])
                                            big_bin.extend(self.bins[x-1][y+1][z-1])
                                            big_bin.extend(self.bins[x][y-1][z+1])
                                            big_bin.extend(self.bins[x][y+1][z+1])
                                            big_bin.extend(self.bins[x-1][y-1][z+1])
                                            big_bin.extend(self.bins[x-1][y+1][z+1])
                                        else:
                                            if z > 0:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x-1][y][z-1])
                                                big_bin.extend(self.bins[x][y-1][z-1])
                                                big_bin.extend(self.bins[x][y+1][z-1])
                                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                            if z < zmax:
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x-1][y][z+1])
                                                big_bin.extend(self.bins[x][y-1][z+1])
                                                big_bin.extend(self.bins[x][y+1][z+1])
                                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                                big_bin.extend(self.bins[x-1][y+1][z+1])
                                    else:
                                        if y > 0:
                                            big_bin.extend(self.bins[x][y-1][z])
                                            big_bin.extend(self.bins[x-1][y-1][z])
                                            if 0 < z < zmax:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x-1][y][z-1])
                                                big_bin.extend(self.bins[x-1][y][z+1])
                                                big_bin.extend(self.bins[x][y-1][z-1])
                                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                                big_bin.extend(self.bins[x][y-1][z+1])
                                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                            else:
                                                if z > 0:
                                                    big_bin.extend(self.bins[x][y][z-1])
                                                    big_bin.extend(self.bins[x-1][y][z-1])
                                                    big_bin.extend(self.bins[x][y-1][z-1])
                                                    big_bin.extend(self.bins[x-1][y-1][z-1])
                                                if z < zmax:
                                                    big_bin.extend(self.bins[x][y][z+1])
                                                    big_bin.extend(self.bins[x-1][y][z+1])
                                                    big_bin.extend(self.bins[x][y-1][z+1])
                                                    big_bin.extend(self.bins[x-1][y-1][z+1])
                                        if y < ymax:
                                            big_bin.extend(self.bins[x][y+1][z])
                                            big_bin.extend(self.bins[x-1][y+1][z])
                                            if 0 < z < zmax:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x-1][y][z-1])
                                                big_bin.extend(self.bins[x-1][y][z+1])
                                                big_bin.extend(self.bins[x][y+1][z-1])
                                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                                big_bin.extend(self.bins[x][y+1][z+1])
                                                big_bin.extend(self.bins[x-1][y+1][z+1])
                                            else:
                                                if z > 0:
                                                    big_bin.extend(self.bins[x][y][z-1])
                                                    big_bin.extend(self.bins[x-1][y][z-1])
                                                    big_bin.extend(self.bins[x][y+1][z-1])
                                                    big_bin.extend(self.bins[x-1][y+1][z-1])
                                                if z < zmax:
                                                    big_bin.extend(self.bins[x][y][z+1])
                                                    big_bin.extend(self.bins[x-1][y][z+1])
                                                    big_bin.extend(self.bins[x][y+1][z+1])
                                                    big_bin.extend(self.bins[x-1][y+1][z+1])
                                if x < xmax:
                                    big_bin.extend(self.bins[x+1][y][z])
                                    if 0 < y < ymax:
                                        big_bin.extend(self.bins[x][y-1][z])
                                        big_bin.extend(self.bins[x][y+1][z])
                                        big_bin.extend(self.bins[x+1][y-1][z])
                                        big_bin.extend(self.bins[x+1][y+1][z])
                                        if 0 < z < zmax:
                                            big_bin.extend(self.bins[x][y][z-1])
                                            big_bin.extend(self.bins[x][y][z+1])
                                            big_bin.extend(self.bins[x+1][y][z-1])
                                            big_bin.extend(self.bins[x+1][y][z+1])
                                            big_bin.extend(self.bins[x][y-1][z-1])
                                            big_bin.extend(self.bins[x][y+1][z-1])
                                            big_bin.extend(self.bins[x+1][y-1][z-1])
                                            big_bin.extend(self.bins[x+1][y+1][z-1])
                                            big_bin.extend(self.bins[x][y-1][z+1])
                                            big_bin.extend(self.bins[x][y+1][z+1])
                                            big_bin.extend(self.bins[x+1][y-1][z+1])
                                            big_bin.extend(self.bins[x+1][y+1][z+1])
                                        else:
                                            if z > 0:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x+1][y][z-1])
                                                big_bin.extend(self.bins[x][y-1][z-1])
                                                big_bin.extend(self.bins[x][y+1][z-1])
                                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                                big_bin.extend(self.bins[x+1][y+1][z-1])
                                            if z < zmax:
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x+1][y][z+1])
                                                big_bin.extend(self.bins[x][y-1][z+1])
                                                big_bin.extend(self.bins[x][y+1][z+1])
                                                big_bin.extend(self.bins[x+1][y-1][z+1])
                                                big_bin.extend(self.bins[x+1][y+1][z+1])
                                    else:
                                        if y > 0:
                                            big_bin.extend(self.bins[x][y-1][z])
                                            big_bin.extend(self.bins[x+1][y-1][z])
                                            if 0 < z < zmax:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x+1][y][z-1])
                                                big_bin.extend(self.bins[x+1][y][z+1])
                                                big_bin.extend(self.bins[x][y-1][z-1])
                                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                                big_bin.extend(self.bins[x][y-1][z+1])
                                                big_bin.extend(self.bins[x+1][y-1][z+1])
                                            else:
                                                if z > 0:
                                                    big_bin.extend(self.bins[x][y][z-1])
                                                    big_bin.extend(self.bins[x+1][y][z-1])
                                                    big_bin.extend(self.bins[x][y-1][z-1])
                                                    big_bin.extend(self.bins[x+1][y-1][z-1])
                                                if z < zmax:
                                                    big_bin.extend(self.bins[x][y][z+1])
                                                    big_bin.extend(self.bins[x+1][y][z+1])
                                                    big_bin.extend(self.bins[x][y-1][z+1])
                                                    big_bin.extend(self.bins[x+1][y-1][z+1])
                                        if y < ymax:
                                            big_bin.extend(self.bins[x][y+1][z])
                                            big_bin.extend(self.bins[x+1][y+1][z])
                                            if 0 < z < zmax:
                                                big_bin.extend(self.bins[x][y][z-1])
                                                big_bin.extend(self.bins[x][y][z+1])
                                                big_bin.extend(self.bins[x+1][y][z-1])
                                                big_bin.extend(self.bins[x+1][y][z+1])
                                                big_bin.extend(self.bins[x][y+1][z-1])
                                                big_bin.extend(self.bins[x+1][y+1][z-1])
                                                big_bin.extend(self.bins[x][y+1][z+1])
                                                big_bin.extend(self.bins[x+1][y+1][z+1])
                                            else:
                                                if z > 0:
                                                    big_bin.extend(self.bins[x][y][z-1])
                                                    big_bin.extend(self.bins[x+1][y][z-1])
                                                    big_bin.extend(self.bins[x][y+1][z-1])
                                                    big_bin.extend(self.bins[x+1][y+1][z-1])
                                                if z < zmax:
                                                    big_bin.extend(self.bins[x][y][z+1])
                                                    big_bin.extend(self.bins[x+1][y][z+1])
                                                    big_bin.extend(self.bins[x][y+1][z+1])
                                                    big_bin.extend(self.bins[x+1][y+1][z+1])
        print(sorted(big_bin))

a = Aaa()
a.aaa()