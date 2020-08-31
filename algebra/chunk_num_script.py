# Python script to chunk numbers into 64-bit hexadecimal numbers:
lst = list("cfca638f1500e327035cdf02acb2744d06e68545f7e64c256ab7ae14297a1a823132b971cdefc65870636cb60d217ff87fa59308c07a8fab8579e02ed3cddca5b093ed79b1c57b5fe3f89c11811c1e214983de300000535e7bc00000000060")
def get(lst, n):
     return ['0x' + ''.join(lst[-n:])] + ["0x" + ''.join(lst[(-i-n):-i]) for i in range(n, len(lst), n)]
[print("{},".format(x)) for x in get(lst, 16)]
