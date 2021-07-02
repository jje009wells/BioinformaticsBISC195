module Assignment07

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       isDNA
#=
include("assignment3code.jl")
include("assignment4code.jl")
include("assignment5code.jl")
include("assignment6code.jl")
=#

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    composition(sequence)

    Counts the number of each type of base
    in a DNA sequence and returns a Dict of those counts for A, C, G, T, and N

    Ex.
    julia> composition("ACCGGGTTTTN")
   Dict{Char,Int64} with 5 entries:
     'A' => 1
     'G' => 3
     'T' => 4
     'N' => 1
     'C' => 2
    
   julia> composition("AAX")
   ERROR: Invalid base, X
"""
function composition(sequence)
    sequence = uppercase(sequence)
    bases = Dict('A' => 0, 'G' => 0, 'T' => 0, 'N' => 0, 'C' => 0)

    for base in sequence
        if base == 'A'
            bases['A'] = bases['A'] + 1
        elseif base == 'C'
            bases['C'] = bases['C'] + 1
        elseif base == 'G'
            bases['G'] = bases['G'] + 1
        elseif base == 'T'
            bases['T'] = bases['T'] + 1
        elseif base == 'N'
            bases['N'] = bases['N'] + 1
        else
            throw(ErrorException("Invalid base $base"))
        end
    end

    return bases
end

"""
    gc_content(seq=)

    Calculates the GC ratio of a DNA sequence.
    The GC ratio is the total number of G and C bases divided by the total length of the sequence.
    For more info about GC content, see here: https://en.wikipedia.org/wiki/GC-content

    Ex.

   julia> gc_content("ATNG")
   0.25
   
   julia> gc_content("ccccggggn")
   0.8888888888888888

"""
function gc_content(seq)
    if !isDNA(seq)
        throw(ErrorException("Invalid sequence $seq"))
    end
    seq = uppercase(seq)

    seqlength = length(seq)

    ## count the number of G's
    gs = count(==('G'), seq)
    ## count the number of C's
    cs = count(==('C'), seq)

    return (gs + cs) / seqlength
end

"""
    complement(base)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `Char`.
If a valid base is not provided, the function throws an error.
"""
function complement(base::Char)
    complements = Dict("A" => 'T',
                       "T" => 'A',
                       "G" => 'C',
                       "C" => 'G')
    
    base = uppercase(string(base))
    
    !(base in keys(complements)) && error("Invalid base $base")
    return complements[base]
end

"""
    complement(base)

    Get the DNA complement of the provided sequence

    Ex.
    julia> complement("ATTN")
    "TAAN"

    julia> complement("ATTAGC")
    "TAATCG"

    Accepts Strings only, returns string
    If a valid base is not provided, the function throws an error.
"""
function complement(seq::String)
    seq = uppercase(seq)
    returnSeq = ""
    for base in seq
        if base == 'N'
            returnSeq = returnSeq * base
        else
            returnSeq = returnSeq * complement(base)
        end
    end

    return returnSeq

end

"""
    reverse_complement(sequence)

Takes a DNA sequence and returns the reverse complement
of that sequence.

Takes lowercase or uppercase sequences,
but always returns uppercase.

Examples
≡≡≡≡≡≡≡≡≡≡
    julia> reverse_complement("AAATTT")
    "AAATTT"

    julia> reverse_complement("GCAT")
    "ATGC"

    julia> rc = reverse_complement("TTGGG");

    julia> println(rc)
    CCCAA
    
    julia> reverse_complement("ATTAGC")
    "GCTAAT"

    julia> reverse_complement("ATN")
    "NAT"
"""
function reverse_complement(sequence)
    reverseseq = complement(sequence)
    reverseseq = reverse(reverseseq)
    return reverseseq
end

"""
 function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing headers,
the other containing the entire sequence as a `String`.

Note: function validates DNA sequences for correctness.

Example
≡≡≡≡≡≡≡≡≡
julia> ex1 = parse_fasta("data/ex1.fasta");

julia> ex1[1]
2-element Array{String,1}:
   "ex1.1 | easy"
   "ex1.2 | multiline"

julia> ex1[2]
2-element Array{String,1}:
   "AATTATAGC"
   "CGCCCCCCAGTCGGATT"

julia> ex2 = parse_fasta("data/ex2.fasta");
ERROR: invalid base H
"""
function parse_fasta(path)
 headersVect = Vector()
 seqVect = Vector()
 toBeCombined = Vector()

 for line in eachline(path)
     #@info "line is " line
     #@info "to be combined is" toBeCombined
     if startswith(line, '>')
         push!(headersVect, line[2:end])
         if size(toBeCombined,1) > 0
             push!(seqVect, join(toBeCombined))
             toBeCombined = Vector()
         end
     elseif !isDNA(line)
        throw(ErrorException("Invalid base $base"))
     elseif isDNA(line)
         push!(toBeCombined, line)

     end
     #@info "to be combined is" toBeCombined
 end
 push!(seqVect, join(toBeCombined))
 
 return (headersVect, seqVect)
end

"""
    isDNA(sequence)

checks the given sequence, and if one of its letters is not a base as stored in the bases array, then returns false.
otherwise, returns true
"""
function isDNA(sequence)
    bases = ['A', 'C', 'G', 'T', 'N']
    sequence = uppercase(sequence)
    for letter in sequence
        if !(letter ∈ bases)
            return false
        end
    end
    return true
end
   # Your code here.
# Don't forget to export your functions!


end # module Assignment07
