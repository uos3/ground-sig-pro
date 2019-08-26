import math


def get_allowed_k_values():
    out = []
    out += list(range(40,512+1,8))
    out += list(range(512,1024+1,16))
    out += list(range(1024,2048+1,32))
    out += list(range(2048,6144+1,64))
    return out

def get_k_index(k):
    allowed = get_allowed_k_values()
    return allowed.index(k)

def get_f1_f2(k_index):
    # from 5.1.3.2.3
    f1_list = [3, 7, 19, 7, 7, 11, 5, 11, 7, 41, 103, 15, 9, 17, 9, 21, 101, 21, 57, 23, 13, 27, 11, 27, 85, 29, 33, 15, 17, 33, 103, 19, 19, 37, 19, 21, 21, 115, 193, 21, 133, 81, 45, 23, 243, 151, 155, \
               25, 51, 47, 91, 29, 29, 247, 29, 89, 91, 157, 55, 31, 17, 35, 227, 65, 19, 37, 41, 39, 185, 43, 21, 155, 79, 139, 23, 217, 25, 17, 127, 25, 239, 17, 137, 215, 29, 15, 147, 29, 59, 65, 55, 31, 17, 171, \
               67, 35, 19, 39, 19, 199, 21, 211, 21, 43, 149, 45, 49, 71, 13, 17, 25, 183, 55, 127, 27, 29, 29, 57, 45, 31, 59, 185, 113, 31, 17, 171, 209, 253, 367, 265, 181, 39, 27, 127, 143, 43, 29, 45, 157, 47, 13, \
               111, 443, 51, 51, 451, 257, 57, 313, 271, 179, 331, 363, 375, 127, 31, 33, 43, 33, 477, 35, 233, 357, 337, 37, 71, 71, 37, 39, 127, 39, 39, 31, 113, 41, 251, 43, 21, 43, 45, 45, 161, 89, 323, 47, 23, 47, 263]
    f2_list = [10, 12, 42, 16, 18, 20, 22, 24, 26, 84, 90, 32, 34, 108, 38, 120, 84, 44, 46, 48, 50, 52, 36, 56, 58, 60, 62, 32, 198, 68, 210, 36, 74, 76, 78, 120, 82, 84, 86, 44, 90, 46, 94, 48, 98, 40, 102, \
               52, 106, 72, 110, 168, 114, 58, 118, 180, 122, 62, 84, 64, 66, 68, 420, 96, 74, 76, 234, 80, 82, 252, 86, 44, 120, 92, 94, 48, 98, 80, 102, 52, 106, 48, 110, 112, 114, 58, 118, 60, 122, 124, 84, 64, 66, 204, \
               140, 72, 74, 76, 78, 240, 82, 252, 86, 88, 60, 92, 846, 48, 28, 80, 102, 104, 954, 96, 110, 112, 114, 116, 354, 120, 610, 124, 420, 64, 66, 136, 420, 216, 444, 456, 468, 80, 164, 504, 172, 88, 300, 92, 188, 96, 28, \
               240, 204, 104, 212, 192, 220, 336, 228, 232, 236, 120, 244, 248, 168, 64, 130, 264, 134, 408, 138, 280, 142, 480, 146, 444, 120, 152, 462, 234, 158, 80, 96, 902, 166, 336, 170, 86, 174, 176, 178, 120, 182, 184, 186, 94, 190, 480]
             
    return (f1_list[k_index], f2_list[k_index])
    

class bcjr_decoder(object):
    state_cnt = 8
    trans_cnt = 16
    
    trans_lut = [
    [1-1,        1-1,        0,          0],
    [2-1,        5-1,        0,          0],
    [3-1,        6-1,        0,          1],
    [4-1,        2-1,        0,          1],
    [5-1,        3-1,        0,          1],
    [6-1,        7-1,        0,          1],
    [7-1,        8-1,        0,          0],
    [8-1,        4-1,        0,          0],
    [1-1,        5-1,        1,          1],
    [2-1,        1-1,        1,          1],
    [3-1,        2-1,        1,          0],
    [4-1,        6-1,        1,          0],
    [5-1,        7-1,        1,          0],
    [6-1,        3-1,        1,          0],
    [7-1,        4-1,        1,          1],
    [8-1,        8-1,        1,          1]]
    
    def maxstar(a, b):
        #return max(a,b) + math.log(1.0+math.exp(-abs(a-b)))
        return max(a,b)

    
    def __init__(self, chan_uncoded_llrs, chan_coded_llrs, uncoded_llrs_term, coded_llrs_term):
        
        self.chan_uncoded_llrs = chan_uncoded_llrs
        self.chan_coded_llrs = chan_coded_llrs
        
        self.K = len(self.chan_uncoded_llrs)
        
        # TODO: termination properly
        self.betas_init = [0]*self.state_cnt
    
    def activate(self, interleaved_llrs):
        K = self.K
        state_cnt = self.state_cnt
        trans_cnt = self.trans_cnt
        gammas = [[0] * trans_cnt for i in range(K)]
        alphas = [[0] * state_cnt for i in range(K)]
        betas = [[0] * state_cnt for i in range(K)]
        tmp1 = [0]*state_cnt
        
        # Add the uncoded LLRs
        uncoded = self.chan_uncoded_llrs.copy()
        for k in range(0,K):
            uncoded[k] += interleaved_llrs[k]
        
        # Calculate the gammas
        for k in range(0,K):
            for t in range(0, trans_cnt):
                if (self.trans_lut[t][2] and self.trans_lut[t][3]):
                    gammas[k][t] = uncoded[k] + self.chan_coded_llrs[k]
                elif self.trans_lut[t][2]:
                    gammas[k][t] = uncoded[k]
                elif self.trans_lut[t][3]:
                    gammas[k][t] = self.chan_coded_llrs[k] 
                #else:
                #    gammas[k][t] = 0

        
        # Do the forwards recursion
        # set the initial value
        alphas[0] = [0] + [-9999999]*(state_cnt-1)
        for k in range(1, K):
            for t in range(0, state_cnt):
                tmp1[self.trans_lut[t][1]] = alphas[k-1][self.trans_lut[t][0]] + gammas[k-1][t]
            for t in range(state_cnt, trans_cnt):
                tmp2 = alphas[k-1][self.trans_lut[t][0]] + gammas[k-1][t]
                alphas[k][self.trans_lut[t][1]] = bcjr_decoder.maxstar(tmp1[self.trans_lut[t][1]],tmp2)

        # Do the backwards recursion
        betas[K-1] = self.betas_init
        for k in range(K-2,-1,-1):
            for t in range(0, state_cnt):
                tmp1[self.trans_lut[t][0]] = betas[k+1][self.trans_lut[t][1]] + gammas[k+1][t]
            for t in range(state_cnt, trans_cnt):
                tmp2 = betas[k+1][self.trans_lut[t][1]] + gammas[k+1][t]
                betas[k][self.trans_lut[t][0]] = bcjr_decoder.maxstar(tmp1[self.trans_lut[t][0]], tmp2)

        # Calculate the deltas
        deltas = gammas.copy()
        for t in range(0, trans_cnt):
            t1 = self.trans_lut[t][0]
            t2 = self.trans_lut[t][1]
            for k in range(0,K):
                deltas[k][t] = alphas[k][t1] + betas[k][t2]
        
        # Extrinsic LLR out
        extrinsic_out = [0]*K
        for k in range(0, K):
            p1 = None
            p0 = None
            for t in range(0, trans_cnt):
                if self.trans_lut[t][2]:
                    if not p1 is None:
                        p1 = bcjr_decoder.maxstar(p1,deltas[k][t])
                    else:
                        p1 = deltas[k][t]
                else:
                    if not p0 is None:
                        p0 = bcjr_decoder.maxstar(p0,deltas[k][t])
                    else:
                        p0 = deltas[k][t]
            extrinsic_out[k] = p1-p0-uncoded[k]

        return extrinsic_out
        
def interleave(a):
    # from 5.1.3.2.3

    # pre-allocate the output
    K = len(a)
    k_index = get_k_index(K)
    b = [0]*K
    
    # initialise the QPP interleaver address generator
    (f1,f2) = get_f1_f2(k_index)
    
    g = f1+f2
    f = 0

    for i in range(0, K):
        b[f] = a[i]
        
        # increment the interleaver address generator
        f = f + g
        g = g + f2
        if g >= K:
            g = g - K
        if f >= K:
            f = f - K
    
    return b

def component_encoder(c):
    state = 0
    k = len(c)
    z = [0]*k
    for i in range(0,k):
        state_nxt = state << 1
        if c[i]:
            state_nxt |= ((state >> 1) & 1) ^ ((state >> 2) & 1) ^ 1;
        else:
            state_nxt |= ((state >> 1) & 1) ^ ((state >> 2) & 1);
        if ((state_nxt & 1) ^ (state & 1) ^ ((state >> 2) & 1)):
            z[i] = 1
        state = state_nxt
    
    # get the termination
    if  state&7 == 0:
        x_term = [0,0,0]
        z_term = [0,0,0]
    elif state&7 == 1:
        x_term = [0,1,1]
        z_term = [1,0,1]
    elif state&7 == 2:
        x_term = [1,1,0]
        z_term = [0,1,0]
    elif state&7 == 3:
        x_term = [1,0,1]
        z_term = [1,1,1]
    elif state&7 == 4:
        x_term = [1,0,0]
        z_term = [1,0,0]
    elif state&7 == 5:
        x_term = [1,1,1]
        z_term = [0,0,1]
    elif state&7 == 6:
        x_term = [0,1,0]
        z_term = [1,1,0]
    else: #state&7 == 7:
        x_term = [0,0,1]
        z_term = [0,1,1]
        
    return (z, x_term, z_term)

def subblock_interleave(d, d_select):
    D = len(d)
    C = 32
    R = math.ceil(D/C)
    N_D = R*C-D  # number of prepended nulls
    y = [-1]*N_D +  d
    P =  [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31]

    v = [0]*(R*C)
    
    if not d_select == 2:
        # Perform interleaving for d_0, d_1
        # Read out column wise, jumping between columns based on the pattern P
        k = 0
        for c in range(0, C):
            col = P[c]
            for r in range(0, R):
                v[k] = y[ r*C + col ]
                k += 1
    else:
        # For d_2, just the equation
        K = R*C
        for k in range(0, R*C):
            pi_k = ( P[math.floor(k/R)] + C*(k%R) + 1 )%K
            v[k] = y[pi_k]
    return v

def pick_k(A):
    K_list = get_allowed_k_values()
    K = -1
    for k in range(0, len(K_list)):
        if A <= K_list[k]:
            K = K_list[k]
            break
    if (K < 0):
        # error
        return None
    return K

def bit_collection(v_0, v_1, v_2):
    # 5.1.4.1.2 bit collection
    w = v_0
    for i in range(0, len(v_1)):
        w += [v_1[i], v_2[i]]
    return w

def remove_nulls(a):
    b = []
    for i in range(0, len(a)):
        if not a[i] is None:
            b += [a[i]]
    return b
    
def encoder_core(c, F):

    k = len(c)
    
    c_interleaved = interleave(c)
    
    (z, x_term, z_term) = component_encoder(c)
    (z_prime, x_prime_term, z_prime_term) = component_encoder(c_interleaved)
    
    # 5.1.3.2.1 locate the bits in the correct place on the three streams
    d_0 = c
    d_1 = z
    d_2 = z_prime
    
    for f in range(0, F):
        d_0[f] = None
        d_1[f] = None
        
    # as in 5.1.3.2.2
    d_0 += [x_term[0], z_term[1], x_prime_term[0], z_prime_term[1]]
    d_1 += [z_term[0], x_term[2], z_prime_term[0], x_prime_term[2]]
    d_2 += [x_term[1], z_term[2], x_prime_term[1], z_prime_term[2]]
    
    return (d_0, d_1, d_2)
    

def codeblock_encoder_chain(a_bits):
    A = len(a_bits)
    K = pick_k(A)
    F = K-A
    c = [0]*F + a_bits
    
    (d_0, d_1, d_2) = encoder_core(c, F)
    
    v_0 = subblock_interleave(d_0, 0)
    v_1 = subblock_interleave(d_1, 1)
    v_2 = subblock_interleave(d_2, 2)
    
    w = bit_collection(v_0, v_1, v_2)
        
    # for ratematching, just remove the NULLs. (Assumes k_0 = 0)
    # The calling function can select a sub-set of bits if it wishes
    # e should have a length of 3K+12-2F
    e = remove_nulls(w)
            
    return e
    
def decoder_core(d_0, d_1, d_2):
    ## TODO: need to input max iterations; needs CRC early stopping
    K = len(d_0)-4
    
    iterations = 10
    
    sys = d_0[0:K]
    parity_upper = d_1[0:K]
    parity_lower = d_2[0:K]
    sys_interleaved = interleave(sys)
    
    interleave_pattern = interleave(list(range(0,K)))
    
    ##TODO: termination
    
    Upper = bcjr_decoder(sys, parity_upper, 0, 0)
    Lower = bcjr_decoder(sys_interleaved, parity_lower, 0, 0)
    
    upper_ap = [0]*K
    lower_ap = [0]*K

    for i in range(0,iterations):
        upper_ex = Upper.activate(upper_ap)
        for k in range(0,K):
            lower_ap[k] = upper_ex[interleave_pattern[k]]
        lower_ex = Lower.activate(lower_ap)
        for k in range(0,K):
            upper_ap[interleave_pattern[k]] = lower_ex[k]
            
    # Now we're done, output the final LLR output
    llrs_out = [0]*K
    for k in range(0,K):
        llrs_out[k] = upper_ap[k] + upper_ex[k] + sys[k]
    return llrs_out
    
    
def codeblock_decoder_chain(e_llrs, A):
    K = pick_k(A)
    F = K-A
    
    # This jumps over bit collection and subblock interleaving
    #  in one step, by using one interleaver pattern
    indicies_d_0 = list(range(0,K+4))
    indicies_d_1 = list(range(10000,10000+K+4))
    indicies_d_2 = list(range(20000,20000+K+4))
    for f in range(0, F):
        indicies_d_0[f] = None
        indicies_d_1[f] = None
    indicies_v_0 = subblock_interleave(indicies_d_0, 0)
    indicies_v_1 = subblock_interleave(indicies_d_1, 1)
    indicies_v_2 = subblock_interleave(indicies_d_2, 2)
    indicies_w = bit_collection(indicies_v_0, indicies_v_1, indicies_v_2)
    indicies_e = remove_nulls(indicies_w)
    
    # Now use the indicies to jump from e to d
    d_0 = [0]*(K+4)
    d_1 = [0]*(K+4)
    d_2 = [0]*(K+4)
    for i in range(0 ,len(e_llrs)):
        index = indicies_e[i]
        if index < 10000:
            d_0[index] = e_llrs[i] # (change to += if repetition is ever supported)
        elif index < 20000:
            d_1[index-10000] = e_llrs[i]
        else:
            d_2[index-20000] = e_llrs[i]
            
    c_hat = decoder_core(d_0, d_1, d_2)
    
    a_hat = c_hat[F:K-1]
    return a_hat