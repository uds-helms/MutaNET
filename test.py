from tools import SubstitutionMatrix

matrix = SubstitutionMatrix('substitution_matrices/PAM10.tsv')

"""
mutations:
    1. missense:    point mutation
    2. frameshift:  completely different starting with mutation
    3. readthrough: later stop
    4. nonsense:    early stop
    
    
score(ref, mut, pos=-1):
    case 1: mut and ref are empty
        return 1.0
    case 2: ref is empty and mut not
        return 0.0
    case 3: mut is empty and ref not
        return 0.0
        
    case 4: ref and mut not same length
        case 4.1: pos is given
            assumption: ref and mut are the same up to the mutation
            => ref_aligned, mut_aligned = align ref[pos:], mut[pos:]
            => ref = ref[:pos] + ref_aligned
            => mut = mut[:pos] + mut_aligned
        case 4.2: pos is not given
            fill up the shorter of the two with gaps (*) at the end
            => have the same length
    
    case 5: ref and mut have the same length
    
=> case 4.1 must return something of unequal length
            > ref[:pos] and mut[:pos] should have the same length
            
            
"""

ref = 'AAAAAAAAA'
mut = 'RRRRR'

print(matrix.score(ref, mut, 6))

