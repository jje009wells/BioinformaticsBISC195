using BioinformaticsBISC195
using Test

@testset "BioinformaticsBISC195" begin
    
@testset "Using Strings" begin
    
    @testset "normalizeDNA" begin
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws Exception normalizeDNA("ZCA")
        @test_throws Exception normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition" begin
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Dict

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)
    end # composition

    @testset "gc_content" begin
        @test gc_content("ANTG") == 1/3
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
    end # gc_content

    @testset "complement" begin
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception complement("ABC")
    end # complement

    @testset "reverse_complement" begin
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception reverse_complement("ABC")
    end # reverse_complement

    @testset "parse_fasta" begin
        #testpath = normpath(joinpath(@__DIR__, "..", "BioinformaticsBISC195/data"))
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        @test_throws Exception parse_fasta(ex2_path)

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

    @testset "fasta_header" begin
        @test fasta_header("M0002 |China|Homo sapiens") == ("M0002", "China", "Homo sapiens")
        @test fasta_header("Another sequence") == ("Another sequence",)
        @test fasta_header("headers| can | have| any number | of | info blocks") == ("headers", "can", "have", "any number", "of", "info blocks")       
    end

    #test isDNA, kmercount, kmerdistance
    @testset "isDNA" begin
        @test isDNA("AcGt") == true
        @test isDNA("Ackk") == false #returns false bc the normalizeDNA is not used to fix it beforehand
        @test isDNA(normalizeDNA("Ackk")) == true # returns true because normalizeDNA converts the 'k's into 'N's, which is a valid ambiguous base
        @test isDNA("zzcGt") == false
    end

    @testset "kmercount" begin
        @test kmercount("ggg", 3) == Dict("GGG" => 1)
        @test kmercount("ATATATATA", 4) == Dict("TATA" => 3, "ATAT" => 3)
        @test_throws Exception kmercount("ACGGGqq", 4)
        @test_throws Exception kmercount("a", 2)
    end

    @testset "kmerdistance" begin
        @test kmerdistance(["AAA", "AGT"], ["AAA", "GGG"]) == 1 - 1/3
        @test kmerdistance(["AAA", "AGT"], ["TTT", "GGG"]) == 1
        @test kmerdistance(["AAA", "ATT"], ["AGT", "AAA", "GTGG"]) == 0.75
        @test kmerdistance(["GGT", "CAT"], ["CAT", "GGT"]) == 0
    end

    @testset "kmercollecting" begin
        ex1 = parse_fasta(normpath(joinpath(@__DIR__, "..", "data/ex1.fasta")))
        #ex1 = parse_fasta(joinpath(@__DIR__, "data", "ex1.fasta"))
        kc1 = kmercollecting(ex1[2], 3)
        # I changed these all to Sets because the order of a set does not matter for comparison
        v1 = Set(["AAT", "ATT", "TAT", "ATA", "TAG", "TTA", "AGC"])
        v2 = Set(["AGT", "CCC", "CGG", "ATT", "GAT", "GCC", "CAG", "GTC", "CGC", "GGA", "CCA", "TCG"])
        @test kc1[1] == v1
        @test kc1[2] == v2

        set2 = ("ACCGGTT", "AAA")
        kc2 = kmercollecting(set2, 3)
        v3 = Set(["GGT", "ACC", "GTT", "CCG", "CGG"])
        v4 = Set(["AAA"])
        @test kc2[1] == v3
        @test kc2[2] == v4
    end

    @testset "shortendates" begin
        asia = parse_fasta(normpath(joinpath(@__DIR__, "../", "COV2analysis", "data", "covgen_asiashort.fasta")))
        #asia = parse_fasta(joinpath(@__DIR__, "data", "covgen_asiashort.fasta"))
        @test shortendates(asia[1]) == ["2019-12", "2020-01", "2020-02", "2020-03"]

        dates1 = ["Should Work | 2031-04", "Should work|2000-12-01| extra", "should not work | 2020", "should work|2001-01"]
        @test shortendates(dates1) == ["2031-04", "2000-12", "2001-01"]

        dates2 = ["should not work| XX-XX", "should work|XXXX-XX", "should work   |  XXXXXXXextra"] # the function does not check validity of dates, so these should still work
        @test shortendates(dates2) == ["XXXX-XX", "XXXXXXX"]
    end

    @testset "monthlycomparisons" begin
        headers1 = ["Should Work | 2019-12", "Should work|2020-01-01| extra", "should not work | 2020", "should work|2020-02", "should work|2020-03-2"]
        seq1 = ["ACCTTGGA", "ACCCTGGA", "ACCCTGG", "ACCCTGC", "AACCTGC"] 
        data1 = [headers1, seq1]
        @test monthlycomparison(data1[2][1], data1, 3) == [maximum(swscorematrix(data1[2][1], data1[2][2])), maximum(swscorematrix(data1[2][1], data1[2][3])), maximum(swscorematrix(data1[2][1], data1[2][4]))]

        headers2 = ["Should Work | 2019-12", "Should work|2020-01-01| extra", "should skip | 2020-01", "should work|2020-02", "should work|2020-03-2"]
        seq2 = ["CCATAGG", "CCATTGG", "CCATGG", "CAATGG", "CAATGT"]
        data2 = [headers2, seq2]
        @test monthlycomparison(data2[2][1], data2, 3) == [maximum(swscorematrix(data2[2][1], data2[2][2])), maximum(swscorematrix(data2[2][1], data2[2][4])), maximum(swscorematrix(data2[2][1], data2[2][5]))]

        # covgen_asiashort.fasta is a version of the covgen_asia file with only 4 data points listed, from month 2019-21 thru 2020-03
        asia = parse_fasta(normpath(joinpath(@__DIR__, "..", "data", "covgen_asiashort.fasta")))
        #asia = parse_fasta(joinpath(@__DIR__, "data", "covgen_asiashort.fasta"))
        
        @test monthlycomparison(asia[2][1], asia, 3) == [maximum(swscorematrix(asia[2][1], asia[2][2])), maximum(swscorematrix(asia[2][1], asia[2][3])), maximum(swscorematrix(asia[2][1], asia[2][4]))]
        @test_throws Exception monthlycomparison(asia[2][1], asia, 4)
        @test_throws Exception monthlycomparison("Isn't real DNA", asia, 3)

        # this file has the same data as covgen_asiashort, but it skips data for Feb 2020, so it should also throw an error
        asiabad = parse_fasta(normpath(joinpath(@__DIR__, "..", "data", "covgen_asiabad.fasta")))
        #asiabad = parse_fasta(joinpath(@__DIR__, "data", "covgen_asiabad.fasta"))
        @test_throws Exception monthlycomparison(asiabad[2][1], asiabad, 3)
    end

end # strings

# @testset "Using BioSequences" begin
    
#     @testset "normalizeDNA" begin
#         @test normalizeDNA("aatgn") == dna"AATGN"
#         @test_throws Exception normalizeDNA("ZCA")
#         @test_throws Exception normalizeDNA(42)
#         c = normalizeDNA('C') 
#         @test c == dna"c"
#         @test c isa LongSequence
#     end #  normalizeDNA

#     @testset "gc_content" begin
#         @test gc_content(dna"ANTG") == 0.25
#         @test gc_content(dna"cccggg") * 100 == 100.0
#         @test gc_content(dna"ATta") == 0.0
#     end #  composition

#     @testset "composition" begin
#         seq = rand(['A','T','G','C','N'], 20) |> join |> LongDNASeq
#         bc = composition(seq)

#         @test bc[DNA_A] == count(==(DNA_A), collect(seq))
#         @test bc[DNA_C] == count(==(DNA_C), collect(seq))
#         @test bc[DNA_G] == count(==(DNA_G), collect(seq))
#         @test bc[DNA_T] == count(==(DNA_T), collect(seq))
#         @test bc[DNA_N] == count(==(DNA_N), collect(seq))
#     end #  gc_content

#     @testset "complement" begin
#         @test complement(dna"ATTAN") == dna"TAATN"
#         @test complement(dna"gcta") == dna"CGAT"
#         @test complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  complement

#     @testset "reverse_complement" begin
#         @test reverse_complement(dna"ATTAN") == dna"NTAAT"
#         @test reverse_complement(dna"gcta") == dna"TAGC"
#         @test reverse_complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  reverse_complement

#     @testset "parse_fasta" begin
#         testpath = normpath(joinpath(@__DIR__, "..", "data"))
#         genomes = joinpath(testpath, "cov2_genomes.fasta")
#         ex1_path = joinpath(testpath, "ex1.fasta")
#         ex2_path = joinpath(testpath, "ex2.fasta")

#         ex1 = parse_fasta(ex1_path)
#         @test ex1 isa Tuple
#         @test all(x-> x isa String, ex1[1])
#         @test all(x-> x isa LongSequence, ex1[2])

#         @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
#         @test ex1[2] == [dna"AATTATAGC", dna"CGCCCCCCAGTCGGATT"]

#         @test_throws Exception parse_fasta(ex2_path)

#         cov2 = parse_fasta(genomes)
#         @test length(cov2[1]) == 8
#         @test length(cov2[2]) == 8
#     end # parse_fasta

# end # BioSequences

end # BioinformaticsBISC195
