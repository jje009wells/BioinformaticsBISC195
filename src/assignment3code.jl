#=
# start of assignment 3 code
"""
    compliment(base)

Get the DNA compliment of the provided base:

    A <-> T
    G <-> C

Accepts `String` or `Char`, but always returns `Char`.
If a valid base is not provided, the function throws an error.

Examples
≡≡≡≡≡≡≡≡≡≡
    
    julia> compliment('A')
    'T'

    julia> compliment("G")
    'C'

    julia> compliment("T")
    'A'

    julia> compliment('C')
    'G'
"""
function compliment(base)
    base = uppercase(base)
    if base == "A" || base == 'A'
        return 'T'
    elseif base == "T" || base == 'T'
        return 'A'
    elseif base == "G" || base == 'G'
        return 'C'
    elseif base == "C" || base == 'C'
        return 'G'
    else throw(ErrorException("Invalid base provided."))
    end
end


"""
    ispurine(base)

A boolean function that returns `true` if the base is a purine (A or G)
and `false` if it is not.
The function only supports bases A, C, G, and T (throws an error for other values).
Accepts `String` or `Char`.

Examples
=========

    julia> ispurine('A')
    true

    julia> ispurine("C")
    false

    julia> if ispurine("G")
               println("It's a purine!")
           else
               println("It's a pyrimidine!")
           end
    It's a purine!

    julia> ispurine('B')
    Error: "Base B not supported")
"""
function ispurine(base)
    base = uppercase(base)
    if base == "A" || base == 'A' || base == "G" || base == 'G'
        return true
    elseif base == "T" || base == 'T' || base == "C" || base == 'C'
        return false
    else throw(ErrorException("Base is not supported."))
    end
end

"""
    ispyrimidine(base)

A boolean function that returns `true` if the base is a pyrimidine (C or T)
and `false` if it is not.
The function only supports bases A, C, G, and T (throws an error for other values).
Accepts `String` or `Char`.

Examples
=========

    julia> ispyrimidine('G')
    false

    julia> ispyrimidine("T")
    true

    julia> if ispyrimidine("G")
               println("It's a pyrimidine!")
           else
               println("It's a purine!")
           end
    It's a purine!

    julia> ispyrimidine('X')
    Error: "Base X not supported"
"""
function ispyrimidine(base)
    return !ispurine(base)
end


"""
    base_type(base)

Determines whether a base is a purine (A or G) or pyrimidine (T or C),
and returns a `String`.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> base_type("G")
    "purine"

    julia> base_type('C')
    "pyrimidine"

    julia> base_type('Z')
    Error: "Base Z not supported"

    julia> x = base_type('A'); println(x)
    purine
"""
function base_type(base)
    base = uppercase(base)
    if ispyrimidine(base)
        return "pyrimidine"
    elseif ispurine(base)
        return "purine"
    else throw(ErrorException("Base is not supported."))
    end
end
# end of assignment 3 code
=#