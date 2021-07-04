#=
# start of assignment 5 code
"""
    complement(base)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `Char`.
If a valid base is not provided, the function throws an error.
"""
function complement(base)
    complements = Dict("A" => 'T',
                       "T" => 'A',
                       "G" => 'C',
                       "C" => 'G')
    
    base = uppercase(string(base))
    
    !(base in keys(complements)) && error("Invalid base $base")
    return complements[base]
end

"""
    isreverse(word1, word2)

    Checks to see if word1 is the reverse of word2
        Returns true if so, false otherwise
"""
function isreverse(word1, word2)
    if length(word1) != length(word2)
        return false
    end
    i = firstindex(word1)
    j = lastindex(word2)
    while j > 0
        
        #@show i j
        if word1[i] != word2[j]
            return false
        end
        i = nextind(word1, i)
        j = prevind(word2, j)
    end
    true
end
"""
    isreversecomplement(seq1, seq2)

Boolean function that checks whether two DNA seqences
are the reverse complement of one another, irrespective of capitalization.
Returns true if yes, false otherwise.

If any invalid bases are encountered,
or if sequences are different length, throws an error.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> isreversecomplement("aaatttcg", "cgaaattt")
    true

    julia> if isreversecomplement("C", "A")
               println("Yes!")
           else
               println("No!")
           end
    No!

    julia> isreversecomplement("TX", "AG")
    Error: Invalid base X

    julia> isreversecomplement("G", "CC")
    Error: Cannot compare sequuences of different length
"""
function isreversecomplement(seq1, seq2)
    ## need to check lengths, then reverse, as reverse if wrong base found then error
    seq2 = uppercase(seq2)
    if length(seq1) != length(seq2)
        throw(ErrorException("Cannot compare sequences of different length"))
    end

    #first find comp of 1, then find out if that is reverse of 2
    seq1comp = ""
    for letter in seq1
        seq1comp = seq1comp * complement(letter)
    end

    return isreverse(seq1comp, seq2)
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
"""
function reverse_complement(sequence)
    ## your code here
    reverseseq = ""
    for letter in sequence
        reverseseq = reverseseq * complement(letter)
    end

    reverseseq = reverse(reverseseq)
end

"""
    isDNA(sequence)

checks the given sequence, and if one of its letters is not a base as stored in the bases array, then returns false.
otherwise, returns true
"""
function isDNA(sequence)
    bases = ['A', 'C', 'G', 'T']
    sequence = uppercase(sequence)
    for letter in sequence
        if !(letter ∈ bases)
            return false
        end
    end
    return true
end
# end of assignment 5 code

"""
   gc_content(sequence)

   Calculates the GC ratio of a DNA sequence.
   The GC ratio is the total number of G and C bases divided by the total length of the sequence.
   For more info about GC content, see here: https://en.wikipedia.org/wiki/GC-content
   
   Example
   ≡≡≡≡≡≡≡≡≡≡
   
       julia> question3("AATG")
       0.25
   
       julia> question3("CCCGG")
       1.0
   
       julia> question3("ATTA")
       0.0
"""
function gc_content(sequence)
    sequence = uppercase(sequence)

    seqlength = length(sequence)

    ## count the number of G's
    gs = count(==('G'), sequence)
    ## count the number of C's
    cs = count(==('C'), sequence)

    return (gs + cs) / seqlength
end
=#