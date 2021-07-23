module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       isDNA,
       kmercount,
       kmerdistance,
       kmercollecting,
       kmercombining,
       monthlycomparison,
       shortendates,
       fasta_header,
       nwscore,
       nwsetupmatrix,
       nwscorematrix,
       swsetupmatrix,
       swscorematrix,
       nwalign,
       swalign
#=
include("assignment3code.jl")
include("assignment4code.jl")
include("assignment5code.jl")
include("assignment6code.jl")
=#

include("AlignmentAlgs.jl")


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
    seq = replace(seq, "Y" => "N")
    seq = replace(seq, "W" => "N")
    seq = replace(seq, "K" => "N")
    seq = replace(seq, "M" => "N")
    seq = replace(seq, "B" => "N")
    seq = replace(seq, "D" => "N")
    seq = replace(seq, "H" => "N")
    seq = replace(seq, "V" => "N")
    seq = replace(seq, "R" => "N")
    seq = replace(seq, "S" => "N")
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
    Ambiguous bases (N) not included in calculations.

    Ex.

   julia> gc_content("ATNG")
   0.25
   
   julia> gc_content("ccccggggn")
   0.8888888888888888

"""
function gc_content(seq)
    seq = normalizeDNA(seq)
    if !isDNA(seq)
        throw(ErrorException("Invalid sequence $seq"))
    end
    seq = normalizeDNA(seq)

    #seqlength = length(seq)

    ## count the number of G's
    gs = count(==('G'), seq)
    ## count the number of C's
    cs = count(==('C'), seq)

    ts = count(==('T'), seq)
    as = count(==('A'), seq)

    return (gs + cs) / (gs + cs + ts + as)
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
 headersVect = Vector{String}()
 seqVect = Vector{String}()
 toBeCombined = Vector{String}()

 for line in eachline(path)
     #@info "line is " line
     #@info "to be combined is" toBeCombined
     if startswith(line, '>')
         push!(headersVect, line[2:end])
         if size(toBeCombined,1) > 0
             push!(seqVect, join(toBeCombined))
             toBeCombined = Vector()
         end
     elseif !isDNA(normalizeDNA(line))
        throw(ErrorException("Invalid base $base"))
     elseif isDNA(normalizeDNA(line))
         push!(toBeCombined, line)

     end
     #@info "to be combined is" toBeCombined
 end
 push!(seqVect, join(toBeCombined))
 
 return (headersVect, seqVect)
end

"""
    isDNA(sequence)

Examples
========

    isDNA("atGG")
    true

    isDNA("ATYYG")
    false

Checks the given sequence, and if one of its letters is not a base as stored in the bases array, then returns false.
otherwise, returns true.
NOTE: does not make use of normalizeDNA to correct the DNA into a uniform format with accepted bases; that will need to be done before isDNA is used
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
        kmer = normalizeDNA(kmer) 
        if !occursin("N", kmer)
            if haskey(kmers, kmer)
            #       add 1 to the value referenced by that kmer
                kmers[kmer] = kmers[kmer] + 1
            #   otherwise
            else
            #       make a new entry in the dictionary with a value of 1
                kmers[kmer] = 1
                #setindex!(kmers, 1, kmer)
            end
        end
        #println(kmers)
    end
    return kmers
end

"""
    kmerdistance(set1, set2)

Takes two kmer sets and calculates the distance between them, using distance metric that is 1 - (length of intersection / length of union).
Or alternatively, ((length of set1 - set2) + (length of set2 - set1)) / length of union.
Note: does not check for validity of kmers

Examples
========

kmerdistance(["AAA", "AGT"], ["AAA", "GGG"])
0.666667 (2/3)

kmerdistance(["AAA", "AGT"], ["TTT", "GGG"])
1

kmerdistance(["AAA", "ATT"], ["AGT", "AAA", "GTGG"])
0.75

"""
function kmerdistance(set1, set2)
    intersectLength = length(intersect(set1, set2))
    unionLength = length(union(set1, set2))
    return 1 - (intersectLength / unionLength)
end


"""
    kmercollecting(set, k)
    
    Takes in a set of sequences and k, then returns a collection of vectors containing all kmers of said length k in the entire set.
    For example, it can take all sequences that were parsed from a certain file, and it returns all kmers of length k that exist in that file.

    Examples
    ========
    ex1 = parse_fasta("C:/Users/Jen/Documents/my_repo/BioinformaticsBISC195/data/ex1.fasta")
    (Any["ex1.1 | easy", "ex1.2 | multiline"], Any["AATTATAGC", "CGCCCCCCAGTCGGATT"])

    kmercollecting(ex1[2], 3)
    2-element Vector{Any}:
    Any["AAT", "ATT", "TAT", "ATA", "TAG", "TTA", "AGC"]
    Any["AGT", "CCC", "CGG", "ATT", "GAT", "GCC", "CAG", "GTC", "CGC", "GGA", "CCA", "TCG"]

    set1 = ("ACCGGTT", "AAA")
    kmercollecting(set1, 3)
    2-element Vector{Any}
    Any["ACC", "CCG", "CGG", "GGT", "GTT"]
    Any["AAA"]

"""
function kmercollecting(set, k)
    kmerdicts = Vector()
	for seq in set
		push!(kmerdicts, collect(keys(kmercount(seq, k))))
	end
	return kmerdicts	
end


#check that the sequences are organized by date
"""
    monthlycomparison(original, allSeq, totalMonths)

    Takes the string of original sequence that all others will be compared against as well as a collection of all the entire parsed fasta file and the total months to be as input.
    Note: at the moment only works if the original string of data is from Dec 2019
    With more time I would modify this to be much more flexible, ie you would be able to put in you own start and end date rather than hard coding the original dates like I did here.
    Note 1: Also depends on the fact that there is at least one datapoint for each month, which there is in the big file of for all sequences in Asia, but is not necessarily always true
    Note 2: the swscorematrix format is very time consuming so for my purposes the totalMonths must be fairly low or the length of the sequences very short in order to reduce comparisons

    Examples
    ========
    covgen_asiashort.fasta is a version of the covgen_asia file with only 4 data points listed, from month 2019-21 thru 2020-03

"""
function monthlycomparison(original, allSeq, totalMonths)
	# stores all of dates for each sequences in its own file for easier traversal later
    if isDNA(normalizeDNA(original)) == false
        throw(ErrorException("Original sequences is not valid DNA"))
    end
    dates = shortendates(allSeq[1]) # I originally tried to put all of the shortendates code inside monthly comparions,
                                    # but I ended up making them two separate functions for easier testing
	#=for headers in allSeq[1]
		#@info headers
		header = fasta_header(headers)
		#@info "header" header[2]
		if length(header[2]) >= 7 #if its shorter than this then it does not have a month listed, so I can't use it and don't need it in the vector
			push!(dates, header[2][1:7])
			@info "header short" header[2][1:7]
		end
		#@info dates
	end=#
	# @info length(dates)
	# @info length(allSeq[1])
	
	# filling out dateIndices to hold ints that represent the indices of the next correct sequences to grab
    # will be used to grab the correct sequences from allSeq to that will be used to score the differences per month
	dateIndices = Vector{Int}()
	# targetMonth and Year are hard coded for the original sequence having a date of Dec 2019
    targetMonth = 1 
	targetYear = 2020
	for i in 1:totalMonths
		targetString = string(targetYear, "-", lpad(targetMonth,2,"0"))
		# @info targetString
		# @info findfirst(targetString .== dates)
        if findfirst(targetString .== dates) === nothing
            throw(ErrorException("Date $targetString does not exist in file; cannot make complete graph."))
        end
		push!(dateIndices, findfirst(targetString .== dates))
		if targetMonth == 12
			targetMonth = 1
			targetYear = targetYear+1
		else
			targetMonth = targetMonth+1
		end
	end
	
	score = Vector{Int}()
	for index in dateIndices
		@info "index" index
		#@info original
		#@info allSeq[2][index]
		#@info "max"  maximum(swscorematrix(original, allSeq[2][index]))
		push!(score, maximum(swscorematrix(original, allSeq[2][index])))
	end
	return score
end	

"""
    shortendates(dates)

    Takes in a collection of parsed header lines and returns a collection of just the dates of these headers, in the same order, using only year-month format.
    Does not save or return any dates that only have a year listed.
    Designed as a helper function for monthlycomparison
    Note: strongly depends on the format of the headers in the fasta file, namely that the date is the second item the fasta header and the date is in the order year-month-day.

    Examples
    ========
    dates1 = ["Should Work | 2031-04", "Should work|2000-12-01| extra", "should not work | 2020", "should work|2001-01"]
    shortendates(dates1)
    3-element Vector{String}:
    ["2031-04", "2000-12", "2001-01"]
"""
function shortendates(parsed_headers)
    shortDates = Vector{String}()
	for headers in parsed_headers
		#@info headers
		header_components = fasta_header(headers)
		#@info "header" header[2]
		if length(header_components[2]) >= 7 #if its shorter than this then it does not have a month listed, so I can't use it and don't need it in the vector
			push!(shortDates, header_components[2][1:7])
			#@info "header short" header[2][1:7]
		end
		#@info dates
	end

    return shortDates
end

"""
 fasta_header(header)

Divides a fasta header into its component parts,
removing any leading or trailing spaces.
NOTE: does not check for or depend on having a '>' to mark as a header

Example
≡≡≡≡≡≡≡≡≡
 julia> fasta_header("M0002 |China|Homo sapiens")
 ("M0002", "China", "Homo sapiens")

 julia> fasta_header("AAATTC")
 Error: Invalid header (headers must start with '>')

 julia> fasta_header("Another sequence")
 ("Another sequence",)

 julia> fasta_header("headers| can | have| any number | of | info blocks")
 ("headers", "can", "have", "any number", "of", "info blocks")
"""
#= function fasta_header(header)
 #startswith(header, '>') || error("Invalid header (headers must start with '>')")
# changed the 2 to a 1 since the parase_fasta cuts the > anyway
# Then you don't need to index at all... `thing[1:end]` is generally === `thing`
    splitVect = split(header, "|")
    returnTuple = (strip(component) for component in splitVect) # this is waaay more efficient than repeatedly destructuring and rebuilding the Tuple   
     # returnTuple = ()
     # for component in splitVect
     #    returnTuple = (returnTuple..., strip(component))
     # end

    return returnTuple
end =#

function fasta_header(header)
   # I ended up not using  returnTuple = (strip(component) for component in splitVect) bc it returned a Base.generator that was hard to deal with,
   # but I don't think this way constantly creates and destroys a tuple over and over again either
   splitVect = split(header, "|")
   return Tuple(map(strip, splitVect))
end
end # module BioinformaticsBISC195
