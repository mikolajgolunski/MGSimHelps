from enum import Enum


class Aaa(Enum):
    a = 1
    b = 2
    c = 3

print(str(Aaa.a.name))
