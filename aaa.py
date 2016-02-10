import numpy as np
import timeit

coord1 = [1.9, 3.8, -7.56]
coord2 = [98.6, -4.09, -0.7]

timeit.timeit("np.linalg.norm(np.array([1.9, 3.8, -7.56]) - np.array([98.6, -4.09, -0.7]))", "import numpy as np; coord1 = [1.9, 3.8, -7.56]; coord2 = [98.6, -4.09, -0.7]", number=10000)

coord1a = np.array(coord1)
coord2a = np.array(coord2)
timeit.timeit("np.linalg.norm(coord1a - coord2a)", "import numpy as np; coord1a = np.array([1.9, 3.8, -7.56]); coord2a = np.array([98.6, -4.09, -0.7])", number=10000)

timeit.timeit("math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)", "import math; coord1 = [1.9, 3.8, -7.56]; coord2 = [98.6, -4.09, -0.7]", number=10000)

timeit.timeit("math.sqrt((coord1a[0] - coord2a[0])**2 + (coord1a[1] - coord2a[1])**2 + (coord1a[2] - coord2a[2])**2)", "import math; import numpy as np; coord1a = np.array([1.9, 3.8, -7.56]); coord2a = np.array([98.6, -4.09, -0.7])", number=10000)
