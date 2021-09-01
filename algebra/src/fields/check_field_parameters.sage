# The following Sage script check the consistency of the following field parameters:
#   1) MODULUS                      must be a prime number   
#   2) MODULUS_BITS                 must be the binary length of MODULUS
#   3) REPR_SHAVE_BITS              must be = length(BigInteger) * 64 - MODULUS_BITS
#   4) R                            must be = 2^(64*BigIntegerLen) mod MODULUS
#   5) R2                           must be = R^2 mod MODULUS
#   6) INV                          must be = (-MODULUS)^(-1) mod 2^64
#   7) TWO_ADICITY                  must be the maximum power of 2 which divides (MODULUS - 1)
#   8) T                            must be = (MODULUS - 1)/2^(TWO_ADICITY)
#   9) T_MINUS_ONE_DIV_TWO          must be = (T - 1)/2
#   10)MODULUS_MINUS_ONE_DIV_TWO    must be = (MODULUS - 1)/2
#   11)GENERATOR                    must have multiplicative order (MODULUS - 1)
#   12)ROOT_OF_UNITY                must have multiplicative order 2^TWO_ADICITY
#
# Open Sage Shell in the corresponding folder and run the command "sage check_field_paramaters sage [filename]".

import re
import sys


#######################################Defining preliminary functions######################################
#Function which pads a hexadecimal representation of a u64 to have length 16
def hex_pad(st):
    l = len(st)
    if l<16:
        for i in range(0,16-l):
            st = "0" + st
    return(st)

##Function which traslates a BigInteger to an integer:
def BigInteger_to_number(string):
    ls_st = re.findall("([0-9a-fA-Fx]+)(?:u64)?\s*,?", string) ## (?:u64)? is ignored thanks to "?:"
    l = len(ls_st)
    #Converting to hexadecimal representation
    for i in range(0,l):
        if not ls_st[i].startswith("0x"):
            ls_st[i] = hex(Integer(ls_st[i]))
            ls_st[i] = ls_st[i].replace("L","")
    output = ""
    for i in range(0,l):
        output = hex_pad(ls_st[i].replace("0x","",1)) + output 
    output = Integer("0x" + output)
    return output

def BigInteger_len(string):
    ls_st = re.findall("([0-9a-fA-Fx]+)(?:u64)?,?", string) ## (?:u64)? is ignored thanks to "?:"
    l = len(ls_st)
    return l


#Functions which returns the two-adicity of a number together with the odd factor
def two_adicity(n):
    x = n
    y = x % 2
    output = 0
    while y == 0:
        x = x//2
        y = x % 2
        output += 1
    return [output,x]

#######################################Reading the values from the file########################
filename = sys.argv[1]

with open(filename) as myfile:
    readfile = myfile.read()

#### Reading the big integers list and extracting names and values
pattern = "const\s+" + "(\w+)" + ":[\w\s:]+=\s*BigInteger\d*\s*\(\s*\[" + "([0-9a-fA-Fxu\s,]+)" + "\]\s*\)"
big_int_ls = re.findall(pattern,readfile)    #####list of couples of the form ('[VARIABLE_NAME]',"[u64],..,[u64]")
big_int_names = [b[0] for b in big_int_ls] 
big_int_values = [BigInteger_to_number(b[1]) for b in big_int_ls]

BigIntegerLen = BigInteger_len(big_int_ls[0][1])

#### Assigning the names to the variables using locals method 
#https://www.pythonpool.com/python-string-to-variable-name/
for s in big_int_names:
    locals()[s] = big_int_values[big_int_names.index(s)]

####Reading the list of usize values
pattern = "const\s+(\w+):\s*(?:u32|u64)\s*=\s*([0-9a-fA-Fx]+)(?:u32|u64)?\s*;"
int_ls = re.findall(pattern,readfile)
int_names = [b[0] for b in int_ls] 
int_values = [Integer(b[1]) for b in int_ls]

#### Assigning the names to the variables using locals method https://www.pythonpool.com/python-string-to-variable-name/
for s in int_names:
    locals()[s] = int_values[int_names.index(s)]



##################################Checking the constistency of the values############################
if is_prime(MODULUS):
    print("Correct. MODULUS is prime")
else:
    print("WARNING! MODULUS IS NOT PRIME!")

#MODULUS_BITS
if MODULUS_BITS == len(bin(MODULUS))-2:
    print("The value of MODULUS_BITS is correct.")
else:
    print("WARNING! THE VALUE OF MODULUS_BITS IS NOT CORRECT!")


#REPR_SHAVE_BITS
if REPR_SHAVE_BITS == 64 * BigIntegerLen - MODULUS_BITS:
        print("The value of REPR_SHAVE_BITS is correct.")
else:
    print("WARNING! THE VALUE OF REPR_SHAVE_BITS IS NOT CORRECT!")


#R and R2
if R == pow(2, 64*BigIntegerLen, MODULUS):
    print("The value of R is correct.")
else:
    print("WARNING! THE VALUE OF R IS NOT CORRECT!")

if R2 == pow(R, 2, MODULUS):
    print("The value of R2 is correct.")
else:
    print("WARNING! THE VALUE OF R2 IS NOT CORRECT!")


#INV
if INV == inverse_mod(-MODULUS, 2**64):
    print("The value of INV is correct.")
else:
    print("WARNING! THE VALUE OF INV IS NOT CORRECT!")


#TWO_ADICITY, T, T_MINUS_ONE_DIV_TWO
s = two_adicity(MODULUS - 1)
if s[0] == TWO_ADICITY:
    print("The value of TWO_ADICITY is correct.")
else:
    print("WARNING! THE VALUE OF TWO_ADICITY IS NOT CORRECT!")

if T == s[1]:
    print("The value of T is correct.")
else:
    print("WARNING! THE VALUE OF T IS NOT CORRECT!")


#T_MINUS_ONE_DIV_TWO
if T_MINUS_ONE_DIV_TWO == (T - 1) // 2:
    print("The value of T_MINUS_ONE_DIV_TWO is correct.")
else:
    print("WARNING! THE VALUE OF T_MINUS_ONE_DIV_TWO IS NOT CORRECT!")

#MODULUS_MINUS_ONE_DIV_TWO
if MODULUS_MINUS_ONE_DIV_TWO == (MODULUS - 1) // 2:
    print("The value of MODULUS_MINUS_ONE_DIV_TWO is correct.")
else:
    print("WARNING! THE VALUE OF MODULUS_MINUS_ONE_DIV_TWO IS NOT CORRECT")

#GENERATOR and ROOT_OF_UNITY
#Converting GENERATOR from Montgomery form
GENERATOR_STANDARD = GENERATOR * inverse_mod(R,MODULUS) % MODULUS

#Checking that GENERATOR has multiplicative order = MODULUS - 1
F = FiniteField(MODULUS)
if F(GENERATOR_STANDARD).multiplicative_order() == MODULUS - 1:
    print("Correct. GENERATOR is a generator of the multiplicative group.")
else:
    print("WARNING! GENERATOR IS NOT A GENERATOR OF THE MULTIPLICATIVE GROUP")

#Writing ROOT_OF_UNITY in standard form from Montgomery
ROOT_OF_UNITY_STANDARD = ROOT_OF_UNITY * inverse_mod(R,MODULUS) % MODULUS

#Checking if the value of ROOT_OF_UNITY a primitive 2^TWO_ADICITY-root
if F(ROOT_OF_UNITY_STANDARD).multiplicative_order() == 2**(TWO_ADICITY):
    print("Correct. ROOT_OF_UNITY is a primitive 2^TWO_ADICITY-root of 1.")
else:
    print("WARNING! ROOT_OF_UNITY IS NOT A PRIMITIVE 2^TWO_ADICITY-ROOT OF 1.")
