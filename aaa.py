try:
    for a in range(20):
        for b in range(87):
            if b > 5:
                raise ValueError
except ValueError:
    print(a)
    print(b)
