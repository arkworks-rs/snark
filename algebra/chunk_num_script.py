# Python script to chunk numbers into 64-bit hexadecimal numbers:
lst = list("9b3af05dd14f6ec619aaf7d34594aabc5ed1347970dec00452217cc900000008508c00000000001")
def get(lst, n):
     return ['0x' + ''.join(lst[-n:])] + ["0x" + ''.join(lst[(-i-n):-i]) for i in range(n, len(lst), n)]
[print("{},".format(x)) for x in get(lst, 16)]
