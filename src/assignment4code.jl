#=
# start of assignment 4 code
"""
    normalizeDNA(sequence)

Normalizes DNA sequences or bases to uppercase `String`s. 
Throws an error if invalid bases are encountered.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> normalizeDNA("aaatg")
    "AAATG"

    julia> normalizeDNA("XTG")
    ERROR: Invalid base X encountered

    julia> normalizeDNA(5)
    ERROR: Invalid base 5 encountered

    julia> normalizeDNA('G')
    "G"

    julia> b = normalizeDNA('G');

    julia> typeof(b)
    String
"""
function normalizeDNA(sequence)
    # 1. make sequence into a String
    
    sequence = string(sequence)
    # 2. make sequence uppercase
    sequence = uppercase(sequence)

    for base in sequence
        if !occursin(base, "ACGT")
            error("Invalid base $base encountered")
        end
    end

    
    return sequence

    # return the normalized sequence
end

"""
    basecomposition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a tuple those counts
in the order A, C, G, T

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> basecomposition("AATCGGG")
    (2, 1, 3, 1)

    julia> basecomposition('C')
    (0, 1, 0, 0)

    julia> A,C,G,T = basecomposition("accgggtttt")
    (1, 2, 3, 4)

    julia> A
    1

    julia> T
    4

    julia> basecomposition("BACCGGGTTTT")
    ERROR: Invalid base B encountered    
"""
function basecomposition(sequence)
    sequence = normalizeDNA(sequence) # make uppercase string, check invalid bases
    a = c = g = t = 0 # sets all 4 variables to `0`

    for base in sequence
        # add 1 to each base as it occurs
        #println(base)
        if base == 'A'
            a = a+1
        elseif base == 'C'
            c = c+1
        elseif base == 'G'
            g = g+1
        elseif base == 'T'
            t = t+1
        end

    end
    
    return a,c,g,t
end

"""
   kmercount(sequence, k)

Finds all kmers in a sequence,
returning a dictionary of those kmers
and the number of times they appear in the sequence.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> kmercount("ggg", 3)
    Dict{Any,Any} with 1 entry:
    "GGG" => 1

    julia> kmercount("ATATATATA", 4)
    Dict{Any,Any} with 2 entries:
    "TATA" => 3
    "ATAT" => 3

    julia> kmercount("ATATATATAx", 4)
    ERROR: Invalid base X encountered

    julia> kmercount("A", 2)
    ERROR: k must be a positive integer less than the length of the sequence
"""
function  kmercount(sequence, k)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Dict() # initialize dictionary
    
    # We're going to loop through the string with numerical index,
    # each time grabbing the bases at position i through i+k-1.
    # What is the last index that we should search?    
    stopindex = length(sequence)-k+1

    for i in 1:stopindex
        kmer = sequence[i:i+k-1] # Change to index the sequence from i to i+k-1
        #println(kmer)
        kmer = normalizeDNA(kmer) 
        #   if this kmer is a key the dictionary
        if haskey(kmers, kmer)
        #       add 1 to the value referenced by that kmer
            kmers[kmer] = kmers[kmer] + 1
        #   otherwise
        else
        #       make a new entry in the dictionary with a value of 1
            kmers[kmer] = 1
            #setindex!(kmers, 1, kmer)
        end
        #println(kmers)

    end
    return kmers
end
# end of assignment 4 code
=#