function swsetupmatrix(s1, s2)
    setupmatrix = zeros(Int, length(s1)+1, length(s2)+1)
    return setupmatrix
end

function swscorematrix(seq1, seq2; match=1, mismatch=-1, gap=-1)
    scoremat = swsetupmatrix(seq1, seq2)
    for i in 2:size(scoremat, 1) # iterate through row indices
        for j in 2:size(scoremat, 2) # iterate through column indices
            above = nwscore(seq1[i-1], nothing; match, mismatch, gap) + scoremat[i-1,j]
            left = nwscore(nothing, seq2[j-1]; match, mismatch, gap) + scoremat[i,j-1]
            diagonal = nwscore(seq1[i-1], seq2[j-1]; match, mismatch, gap) + scoremat[i-1,j-1]
            scoremat[i,j] = max(above, left, diagonal)
            # or `scoremat[i,j] = max(above, left, diagonal, 0)`

            if scoremat[i,j] < 0 
                scoremat[i,j] = 0
            end
        end
    end
    return scoremat
end

"""
    swalign(seq1::String, seq2::String; match =1, mismatch=-1, gap=-1)

    This function takes the results of a matrix scored using smith waterman properties,
    then returns the proper alignment of sequences based on the matrix as a tuple of two strings.

    Ex.)
    julia> swalign("AAACCCGGG","TTTCCCAAA")
    ("CCC", "CCC")

    julia> swalign("TGTTACGG", "GGTTGACTA")
    ("GTT-AC", "GTTGAC")
"""
function swalign(seq1::String, seq2::String; match =1, mismatch=-1, gap=-1)
    if !isDNA(seq1) || !isDNA(seq2) || length(seq1) == 0 || length(seq2) == 0
        throw(ErrorException("Invalid DNA sequence provided"))
    end
    
    scoremat = swscorematrix(seq1, seq2; match=match, mismatch=mismatch, gap=gap)
    aligned1 = ""
    aligned2 = ""
    
    maxScore = findmax(scoremat)

    i = maxScore[2][1]
    j = maxScore[2][2]
    #@info i
    #@info j

    while i > 1 && j > 1
            if scoremat[i,j] == 0
                break
            #check scores from left, right, diag, and if any match current cell, then add that char to the alignments
            elseif scoremat[i,j] == nwscore(seq1[i-1], nothing; match, mismatch, gap) + scoremat[i-1,j]
                #@info "max is from above"
                aligned1 = aligned1 * seq1[i-1]
                aligned2 = aligned2 * "-"
                #update the traversal indices
                i = i-1
            elseif scoremat[i,j] == nwscore(nothing, seq2[j-1]; match, mismatch, gap) + scoremat[i,j-1]
                #@info "max is from left"
                aligned1 = aligned1 * "-"
                aligned2 = aligned2 * seq2[j-1]
                #update the traversal indices
                j = j-1
            else
                #@info "max is from diag"
                aligned1 = aligned1 * seq1[i-1]
                aligned2 = aligned2 * seq2[j-1]
                #update the traversal indices
                i = i-1
                j = j-1
            end

    end
    aligned1 = reverse(aligned1)
    aligned2 = reverse(aligned2)
    return (aligned1, aligned2)
end


# Labs 3 & 4

"""
    nwscore(base1, base2; mismatch = -1, match = 1, gap = -1 )

Your documentation here.

This function compares nucleotides for genomic sequences
to see if they are matching=1, mismatching=-1, or gaps=-1.

Ex.)
```
julia> nwscore('A','G')
-1
julia> nwscore('A','A')
1
julia> nwscore('A',nothing)
-1
julia> nwscore(nothing,nothing)
ERROR: ArgumentError: Score for two gaps is not defined
Stacktrace:
[1] nwscore(#unused#::Nothing, #unused#::Nothing)
    @ Main ~/Documents/my_repo/BISC195Labs/src/needleman_wunch.jl:25
[2] top-level scope
    @ REPL[81]:1
```

"""
function nwscore(base1::Char, base2::Char; match = 1, mismatch = -1, gap = -1)
    if base1 == base2
        return match
    else base1 != base2
        return mismatch
    end
end

function nwscore(base::Char, ::Nothing; match = 1, mismatch = -1, gap = -1)
    return gap
end

function nwscore(::Nothing, base::Char; match = 1, mismatch = -1, gap = -1)
    nwscore(base, nothing; gap=gap)
end

function nwscore(::Nothing, ::Nothing)
    throw(ArgumentError("Score for two gaps is not defined"))
end

# function nwalign(seq1::String, seq2::String; match =1, mismatch=-1, gap=-1)
#     #input sequence 1 (seq1)
#     #input sequence 2 (seq2)
#     #turn each sequence into an array/vector of characters
#     #compare each character based on position (see grid for rules) using nwscore 
#     #make sure there is error if 2 nothings
#     #determine best alignment depending on score
#         #if same score, print both as touple
#     #return alignments (seq1, seq2) as an array
# end 

function nwsetupmatrix(s1, s2; gap=-1)
    setupmatrix = zeros(Int, length(s1)+1, length(s2)+1)
    for j in 1:length(s2)+1
        setupmatrix[1, j] = (j-1)*gap
    end
    for i in 1:length(s1)+1
        setupmatrix[i,1] = (i-1)*gap
    end
    return setupmatrix
end

function nwscorematrix(seq1, seq2; match=1, mismatch=-1, gap=-1)
    scoremat = nwsetupmatrix(seq1, seq2; gap=gap)
    for i in 2:size(scoremat, 1) # iterate through row indices
        for j in 2:size(scoremat, 2) # iterate through column indices
            above = nwscore(seq1[i-1], nothing; match, mismatch, gap) + scoremat[i-1,j]
            left = nwscore(nothing, seq2[j-1]; match, mismatch, gap) + scoremat[i,j-1]
            diagonal = nwscore(seq1[i-1], seq2[j-1]; match, mismatch, gap) + scoremat[i-1,j-1]
            scoremat[i,j] = max(above, left, diagonal)
        end
    end
    return scoremat
end

#need to write a function that finds and returns the best alignment
"""
    nwalign(seq1::String, seq2::String; match =1, mismatch=-1, gap=-1)

    This function takes the results of a matrix scored using needleman_wunsch properties,
    then returns the proper alignment of sequences based on the matrix as a tuple of two strings.

    Ex.)
    julia> nwalign("AGGT", "ACGAT")
    ("AGG-T", "ACGAT")

    julia> nwalign("ACAT", "AGAT")
    ("ACAT", "AGAT")

    julia> nwalign("AGGTC", "ACGT")
    ("AGGTC", "ACGT-")

"""
function nwalign(seq1::String, seq2::String; match =1, mismatch=-1, gap=-1)
    if !isDNA(seq1) || !isDNA(seq2) || length(seq1) == 0 || length(seq2) == 0
        throw(ErrorException("Invalid DNA sequence provided"))
    end
    
    scoremat = nwscorematrix(seq1, seq2; match=match, mismatch=mismatch, gap=gap)
    aligned1 = ""
    aligned2 = ""

    i = size(scoremat,1)
    j = size(scoremat,2)

    while i > 1 && j > 1
            #check scores from left, right, diag, and if any match current cell, then add that char to the alignments
            if scoremat[i,j] == nwscore(seq1[i-1], nothing; match, mismatch, gap) + scoremat[i-1,j]
                #@info "max is from above"
                aligned1 = aligned1 * seq1[i-1]
                aligned2 = aligned2 * "-"
                #update the traversal indices
                i = i-1
            elseif scoremat[i,j] == nwscore(nothing, seq2[j-1]; match, mismatch, gap) + scoremat[i,j-1]
                #@info "max is from left"
                aligned1 = aligned1 * "-"
                aligned2 = aligned2 * seq2[j-1]
                #update the traversal indices
                j = j-1
            else
                #@info "max is from diag"
                aligned1 = aligned1 * seq1[i-1]
                aligned2 = aligned2 * seq2[j-1]
                #update the traversal indices
                i = i-1
                j = j-1
            end

    end
    aligned1 = reverse(aligned1)
    aligned2 = reverse(aligned2)
    return (aligned1, aligned2)
end
