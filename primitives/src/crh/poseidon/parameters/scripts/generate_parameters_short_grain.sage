# This script is based on the corresponding script from https://extgit.iaik.tugraz.at/krypto/hadeshash/-/commit/7ecf9a7d4f37e777ea27e4c4d379443151270563
# We changed the generation of MDS matrices to support short Montgomery multiplication (with respect to the Montgomery constant R2 = 2^64), 
# by sampling random (not Cauchy) MDS matrices with 64 bit short Montgomery representation. 
# For field sizes of 256 bit and beyond the short Montgomery optimization yields better performance than the optimization described in the Appendix B
# of the [Poseidon paper](https://eprint.iacr.org/2019/458). However, regardless of its specific form such short Montgomery MDS matrix should withstand 
# all the attacks discussed in the paper. (needs to be confirmed!)

# Remark: This script contains functionality for GF(2^n), but currently works only over GF(p)! A few small adaptations are needed for GF(2^n).

from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list

# GF(p)= MTN4-753 Fr, x^(-1), N = 2259, n = 753, t = 3, R_F = 8, R_P = 63:
# sage generate_parameters_short_grain.sage 1 1 753 3 8 63 0x1c4c62d92c41110229022eee2cdadb7f997505b8fafed5eb7e8f96c97d87307fdb925e8a0ed8d99d124d9a15af79db26c5c28c859a99b3eebca9429212636b9dff97634993aa4d6c381bc3f0057974ea099170fa13a4fd90776e240000001 > ParametersMNT4Fr.log
# GF(p)= MTN6-753 Fr, x^(-1), N = 2259, n = 753, t = 3, R_F = 8, R_P = 63:
# sage generate_parameters_short_grain.sage 1 1 753 3 8 63 0x1c4c62d92c41110229022eee2cdadb7f997505b8fafed5eb7e8f96c97d87307fdb925e8a0ed8d99d124d9a15af79db117e776f218059db80f0da5cb537e38685acce9767254a4638810719ac425f0e39d54522cdd119f5e9063de245e8001 > ParametersMNT6Fr.log

if len(sys.argv) < 7:
    print("Usage: <script> <field> <s_box> <field_size> <num_cells> <R_F> <R_P> (<prime_number_hex>)")
    print("field = 1 for GF(p)")
    print("s_box = 0 for x^alpha, s_box = 1 for x^(-1)")
    exit()

# Parameters
FIELD = int(sys.argv[1]) # 0 .. GF(2^n), 1 .. GF(p)
SBOX = int(sys.argv[2]) # 0 .. x^alpha, 1 .. x^(-1)
FIELD_SIZE = int(sys.argv[3]) # n
NUM_CELLS = int(sys.argv[4]) # t
R_F_FIXED = int(sys.argv[5])
R_P_FIXED = int(sys.argv[6])

INIT_SEQUENCE = []

PRIME_NUMBER = 0
if FIELD == 1 and len(sys.argv) != 8:
    print("Please specify a prime number (in hex format)!")
    exit()
elif FIELD == 1 and len(sys.argv) == 8:
    PRIME_NUMBER = int(sys.argv[7], 16) # e.g. 0xa7, 0xFFFFFFFFFFFFFEFF, 0xa1a42c3efd6dbfe08daa6041b36322ef

F = GF(PRIME_NUMBER)

def grain_sr_generator():
    bit_sequence = INIT_SEQUENCE
    for _ in range(0, 160):
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        
    while True:
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        while new_bit == 0:
            new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
            new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        yield new_bit
grain_gen = grain_sr_generator()
        
def grain_random_bits(num_bits):
    random_bits = [next(grain_gen) for i in range(0, num_bits)]
    # random_bits.reverse() ## Remove comment to start from least significant bit
    random_int = int("".join(str(i) for i in random_bits), 2)
    return random_int

def init_generator(field, sbox, n, t, R_F, R_P):
    # Generate initial sequence based on parameters
    bit_list_field = [_ for _ in (bin(FIELD)[2:].zfill(2))]
    bit_list_sbox = [_ for _ in (bin(SBOX)[2:].zfill(4))]
    bit_list_n = [_ for _ in (bin(FIELD_SIZE)[2:].zfill(12))]
    bit_list_t = [_ for _ in (bin(NUM_CELLS)[2:].zfill(12))]
    bit_list_R_F = [_ for _ in (bin(R_F)[2:].zfill(10))]
    bit_list_R_P = [_ for _ in (bin(R_P)[2:].zfill(10))]
    bit_list_1 = [1] * 30
    global INIT_SEQUENCE
    INIT_SEQUENCE = bit_list_field + bit_list_sbox + bit_list_n + bit_list_t + bit_list_R_F + bit_list_R_P + bit_list_1
    INIT_SEQUENCE = [int(_) for _ in INIT_SEQUENCE]

def generate_constants(field, n, t, R_F, R_P, prime_number):
    round_constants = []
    num_constants = (R_F + R_P) * t

    if field == 0:
        for i in range(0, num_constants):
            random_int = grain_random_bits(n)
            round_constants.append(random_int)
    elif field == 1:
        for i in range(0, num_constants):
            random_int = grain_random_bits(n)
            while random_int >= prime_number:
                # print("[Info] Round constant is not in prime field! Taking next one.")
                random_int = grain_random_bits(n)
            round_constants.append(random_int)
    return round_constants

def print_round_constants(round_constants, n, field):
    print("Number of round constants:", len(round_constants))

    if field == 0:
        print("Round constants for GF(2^n):")
    elif field == 1:
        print("Round constants for GF(p):")
    hex_length = int(ceil(float(n) / 4)) + 2 # +2 for "0x"
    print(["{0:#0{1}x}".format(entry, hex_length) for entry in round_constants])

# We sample random MDS matrices with short Montgomery representation (with respect
# to the short Montgomery radix R2 = 2^64).
# This approach differs from the recommendation of the Poseidon authors. 
def create_short_mds_p(t):
    while True:
        flag = True
        M = F(2)^-64 * matrix(F,[[F(grain_random_bits(64)) for _ in range(0, t)] for _ in range(0,t)])

	# MDS criterion: all minors of order 1,...,t must be non-zero
	for i in range(1,t+1): 
		if 0 in M.minors(i): 
			flag = False		

        if flag == False:
            continue
        return M

def generate_vectorspace(round_num, M, M_round, NUM_CELLS):
    t = NUM_CELLS
    s = 1
    V = VectorSpace(F, t)
    if round_num == 0:
        return V
    elif round_num == 1:
        return V.subspace(V.basis()[s:])
    else:
        mat_temp = matrix(F)
        for i in range(0, round_num-1):
            add_rows = []
            for j in range(0, s):
                add_rows.append(M_round[i].rows()[j][s:])
            mat_temp = matrix(mat_temp.rows() + add_rows)
        r_k = mat_temp.right_kernel()
        extended_basis_vectors = []
        for vec in r_k.basis():
            extended_basis_vectors.append(vector([0]*s + list(vec)))
        S = V.subspace(extended_basis_vectors)

        return S

def subspace_times_matrix(subspace, M, NUM_CELLS):
    t = NUM_CELLS
    V = VectorSpace(F, t)
    subspace_basis = subspace.basis()
    new_basis = []
    for vec in subspace_basis:
        new_basis.append(M * vec)
    new_subspace = V.subspace(new_basis)
    return new_subspace

# Returns True if the matrix is considered secure, False otherwise
def algorithm_1(M, NUM_CELLS):
    t = NUM_CELLS
    s = 1
    r = floor((t - s) / float(s))

    # Generate round matrices
    M_round = []
    for j in range(0, t+1):
        M_round.append(M^(j+1))

    for i in range(1, r+1):
        mat_test = M^i
        entry = mat_test[0, 0]
        mat_target = matrix.circulant(vector([entry] + ([F(0)] * (t-1))))

        if (mat_test - mat_target) == matrix.circulant(vector([F(0)] * (t))):
            return [False, 1]

        S = generate_vectorspace(i, M, M_round, t)
        V = VectorSpace(F, t)

        basis_vectors= []
        for eigenspace in mat_test.eigenspaces_right(format='galois'):
            if (eigenspace[0] not in F):
                continue
            vector_subspace = eigenspace[1]
            intersection = S.intersection(vector_subspace)
            basis_vectors += intersection.basis()
        IS = V.subspace(basis_vectors)

        if IS.dimension() >= 1 and IS != V:
            return [False, 2]
        for j in range(1, i+1):
            S_mat_mul = subspace_times_matrix(S, M^j, t)
            if S == S_mat_mul:
                print("S.basis():\n", S.basis())
                return [False, 3]

    return [True, 0]

# Returns True if the matrix is considered secure, False otherwise
def algorithm_2(M, NUM_CELLS):
    t = NUM_CELLS
    s = 1

    V = VectorSpace(F, t)
    trail = [None, None]
    test_next = False
    I = range(0, s)
    I_powerset = list(sage.misc.misc.powerset(I))[1:]
    for I_s in I_powerset:
        test_next = False
        new_basis = []
        for l in I_s:
            new_basis.append(V.basis()[l])
        IS = V.subspace(new_basis)
        for i in range(s, t):
            new_basis.append(V.basis()[i])
        full_iota_space = V.subspace(new_basis)
        for l in I_s:
            v = V.basis()[l]
            while True:
                delta = IS.dimension()
                v = M * v
                IS = V.subspace(IS.basis() + [v])
                if IS.dimension() == t or IS.intersection(full_iota_space) != IS:
                    test_next = True
                    break
                if IS.dimension() <= delta:
                    break
            if test_next == True:
                break
        if test_next == True:
            continue
        return [False, [IS, I_s]]

    return [True, None]

# Returns True if the matrix is considered secure, False otherwise
def algorithm_3(M, NUM_CELLS):
    t = NUM_CELLS
    s = 1

    V = VectorSpace(F, t)

    l = 4*t
    for r in range(2, l+1):
        next_r = False
        res_alg_2 = algorithm_2(M^r, t)
        if res_alg_2[0] == False:
            return [False, None]

        # if res_alg_2[1] == None:
        #     continue
        # IS = res_alg_2[1][0]
        # I_s = res_alg_2[1][1]
        # for j in range(1, r):
        #     IS = subspace_times_matrix(IS, M, t)
        #     I_j = []
        #     for i in range(0, s):
        #         new_basis = []
        #         for k in range(0, t):
        #             if k != i:
        #                 new_basis.append(V.basis()[k])
        #         iota_space = V.subspace(new_basis)
        #         if IS.intersection(iota_space) != iota_space:
        #             single_iota_space = V.subspace([V.basis()[i]])
        #             if IS.intersection(single_iota_space) == single_iota_space:
        #                 I_j.append(i)
        #             else:
        #                 next_r = True
        #                 break
        #     if next_r == True:
        #         break
        # if next_r == True:
        #     continue
        # return [False, [IS, I_j, r]]
    
    return [True, None]

def generate_matrix(FIELD, FIELD_SIZE, NUM_CELLS):
    if FIELD == 0:
        print("Matrix generation not implemented for GF(2^n).")
        exit(1)
    elif FIELD == 1:
        mds_matrix = create_short_mds_p(NUM_CELLS)
        result_1 = algorithm_1(mds_matrix, NUM_CELLS)
        result_2 = algorithm_2(mds_matrix, NUM_CELLS)
        result_3 = algorithm_3(mds_matrix, NUM_CELLS)
        while result_1[0] == False or result_2[0] == False or result_3[0] == False:
            mds_matrix = create_short_mds_p(NUM_CELLS)
            result_1 = algorithm_1(mds_matrix, NUM_CELLS)
            result_2 = algorithm_2(mds_matrix, NUM_CELLS)
            result_3 = algorithm_3(mds_matrix, NUM_CELLS)
        return mds_matrix

def print_linear_layer(M, n, t):
    print("n:", n)
    print("t:", t)
    print("N:", (n * t))
    print("Result Algorithm 1:\n", algorithm_1(M, NUM_CELLS))
    print("Result Algorithm 2:\n", algorithm_2(M, NUM_CELLS))
    print("Result Algorithm 3:\n", algorithm_3(M, NUM_CELLS))
    hex_length = int(ceil(float(n) / 4)) + 2 # +2 for "0x"
    print("Prime number:", "0x" + hex(PRIME_NUMBER))
    matrix_string = "["
    for i in range(0, t):
        matrix_string += str(["{0:#0{1}x}".format(int(entry), hex_length) for entry in M[i]])
        if i < (t-1):
            matrix_string += ","
    matrix_string += "]"
    print("MDS matrix:\n", matrix_string)

# Init
init_generator(FIELD, SBOX, FIELD_SIZE, NUM_CELLS, R_F_FIXED, R_P_FIXED)

# Round constants
round_constants = generate_constants(FIELD, FIELD_SIZE, NUM_CELLS, R_F_FIXED, R_P_FIXED, PRIME_NUMBER)
print_round_constants(round_constants, FIELD_SIZE, FIELD)

# Matrix
linear_layer = generate_matrix(FIELD, FIELD_SIZE, NUM_CELLS)
print_linear_layer(linear_layer, FIELD_SIZE, NUM_CELLS)
