# The following Sage script check the consistency of the following curves parameters:
#   1) P=(GENERATOR_X,GENERATOR_Y)  must belongs to the curve of equation E: y^2 = x^3 + Ax + B   
#   2) P                                must have order equal to the MODULUS of the scalar field
#   3) COFACTOR                         must be equal to Order(E)/Order(P)
#   4) COFACTOR_INV                     must be the inverse of COFACTOR in the scalar Field
# Open Sage Shell in the corresponding folder and run the command 
#       "sage check_curve_paramaters sage [file_path_curve] [file_path_basefield] [file_path_scalarfield]".

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
    ls_st = re.findall("([0-9a-fA-Fx]+)(?:u64)?,?", string) ## (?:u64)? is ignored thanks to "?:"
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


#######################################Reading the values from the file containing the curve parameters ########################
filename = sys.argv[1]

with open(filename) as myfile:
    readfile = myfile.read()

##Reading the Base Field and Scalar Field names.
pattern = "type\s*BaseField\s*=\s*(\w+)\s*;"
base_field_name = re.findall(pattern, readfile)[0]
pattern = "type\s*ScalarField\s*=\s*(\w+)\s*;"
scalar_field_name = re.findall(pattern, readfile)[0]
fn = "(?:" + base_field_name + "|" + scalar_field_name + ")" #fn = field name = "(:?Fr|Fq)". Useful declaration for the pattern

#### Reading the big integers list and extracting names and values
pattern = "const\s+(\w+):\s*" + fn + "\s*=\s*field_new!\(\s*" + fn + "\s*,\s*BigInteger\d*\s*\(\s*\[" + "([0-9a-fA-Fxu\s,]+)\s*" + "\]\s*\)"
big_int_ls = re.findall(pattern,readfile)    #####list of couples of the form ('[VARIABLE_NAME]',"[u64],..,[u64]")

big_int_names = [b[0] for b in big_int_ls] 
big_int_values = [BigInteger_to_number(b[1]) for b in big_int_ls]

BigIntegerLen = BigInteger_len(big_int_ls[0][1])

#### Assigning the names to the variables using locals method 
#https://www.pythonpool.com/python-string-to-variable-name/
for s in big_int_names:
    if  "GENERATOR_X" in s: ### Sometimes these variables aren't just defined GENERATOR_X/Y BUT G*_GENERATOR_X/Y
        GENERATOR_X = big_int_values[big_int_names.index(s)]
    elif "GENERATOR_Y" in s:
        GENERATOR_Y = big_int_values[big_int_names.index(s)]
    else:
        locals()[s] = big_int_values[big_int_names.index(s)]

####Reading the value of COFACTOR
pattern =  "const\s+COFACTOR:\s*&'static\s*\[u64\]\s*=\s*&\[([0-9a-fA-Fxu\s,]+)\]\s*;"
COFACTOR = BigInteger_to_number(re.findall(pattern,readfile)[0])

#######################################Reading the values from the file containing the Base Field parameters########################
filename = sys.argv[2]

with open(filename) as myfile:
    readfile = myfile.read()

#### Reading the big integers list and extracting names and values
pattern = "const\s+" + "(\w+)" + ":[\w\s]+=\s*BigInteger\d*\s*\(\s*\[" + "([0-9a-fA-Fxu\s,]+)" + "\]\s*\)"
big_int_ls = re.findall(pattern,readfile)    #####list of couples of the form ('[VARIABLE_NAME]',"[u64],..,[u64]")
big_int_names = [b[0] for b in big_int_ls] 
big_int_values = [BigInteger_to_number(b[1]) for b in big_int_ls]


#### Assigning the names to the variables using locals method 
#https://www.pythonpool.com/python-string-to-variable-name/
for s in big_int_names:
    locals()["BASE_FIELD_" + s] = big_int_values[big_int_names.index(s)]



#######################################Reading the values from the file containing the Scalar Field parameters########################
filename = sys.argv[3]

with open(filename) as myfile:
    readfile = myfile.read()

#### Reading the big integers list and extracting names and values
pattern = "const\s+" + "(\w+)" + ":\s*BigInteger\s*=\s*BigInteger\d*\s*\(\s*\[" + "([0-9a-fA-Fxu\s,]+)" + "\]\s*\)"
big_int_ls = re.findall(pattern,readfile)    #####list of couples of the form ('[VARIABLE_NAME]',"[u64],..,[u64]")
big_int_names = [b[0] for b in big_int_ls] 
big_int_values = [BigInteger_to_number(b[1]) for b in big_int_ls]


#### Assigning the names to the variables using locals method 
#https://www.pythonpool.com/python-string-to-variable-name/
for s in big_int_names:
    locals()["SCALAR_FIELD_" + s] = big_int_values[big_int_names.index(s)]

################################Checking the constistency of the values#############################
Fq = FiniteField(BASE_FIELD_MODULUS)
#Moving from Montgomery to standard form the value of COEFF_A and COEFF_B
A = Fq(COEFF_A) * Fq(BASE_FIELD_R)**(-1)
B = Fq(COEFF_B) * Fq(BASE_FIELD_R)**(-1)
E = EllipticCurve([A,B])
#Moving from Montgomery to standard form the value of COEFF_A and COEFF_B
X = Fq(GENERATOR_X) * Fq(BASE_FIELD_R)**(-1)
Y = Fq(GENERATOR_Y) * Fq(BASE_FIELD_R)**(-1)
if (X,Y) in E:
    print("Correct. (X,Y) is a point of the curve y^2 = x^3 + Ax + B.")
else:
    print("WARNING! (X,Y) IS NOT A POINT OF THE CURVE y^2 = x^3 + Ax + B.")
P = E([X,Y])
###Checking the correctness of the generator
if 0 * P == SCALAR_FIELD_MODULUS * P:
    print("Correct. The MODULUS of the scalar field is equal to the order the generator.")
else:
    print("WARNING! THE MODULUS OF THE SCALAR FIELD IS NOT EQUAL TO THE ORDER OF THE GENERATOR!")
#### Checking the correctness of COFACTOR
if SCALAR_FIELD_MODULUS**2 > 16 * BASE_FIELD_MODULUS:
    # We use the Hasse inequality
    if (BASE_FIELD_MODULUS + 1 - COFACTOR * SCALAR_FIELD_MODULUS)**2 < 4 * BASE_FIELD_MODULUS: 
        print("The value of COFACTOR is correct.")
    else:
        print("WARNING! THE VALUE OF COFACTOR IS NOT CORRECT!")
else:
    # otherwise we use an expensive point counting
    print("Counting the number of points might take long.")
    if COFACTOR * SCALAR_FIELD_MODULUS == E.order():
        print("The value of COFACTOR is correct.")
    else:
        print("WARNING! THE VALUE OF COFACTOR IS NOT CORRECT!")    
###### Checking the correctness of COFACTOR_INV
if Fq(COFACTOR) * Fq(COFACTOR_INV) == Fq(SCALAR_FIELD_R):
    print("Correct. COFACTOR_INV is the inverse of COFACTOR in the the scalar field.")
else:
    print("WARNING! COFACTOR_INV IS NOT THE INVERSE OF COFACTOR IN THE SCALAR FIELD!")