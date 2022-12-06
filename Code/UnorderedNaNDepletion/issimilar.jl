using BioSequences


"Return `true` if a change from `AAᵦ` to `AAᵧ` is considered a similar change according to `aagroups`' classification."
function issimilar(
    AAᵦ::AminoAcid, AAᵧ::AminoAcid,
    aagroups::Dict{AminoAcid,String}
)
    # the two AAs are different
    if AAᵦ ≠ AAᵧ
        # if
        # 1) stop codon isn't classified in any AA group (as expected)
        # 2) one of these AAs is currently a stop codon
        # than these two AAs are considered dissimilar
        AA_Term ∉ keys(aagroups) && (AAᵦ == AA_Term || AAᵧ == AA_Term) && return false
        # return true if, altough the two AAs are different, they still belong in the same group
        return aagroups[AAᵦ] == aagroups[AAᵧ]
        # the two AAs are actually the same -> they must be similar
    else
        return true
    end
end


"Return `true` if a change from `AAᵦ` to `AAᵧ` is considered a similar change according to `substitutionmatrix`."
function issimilar(
    AAᵦ::AminoAcid, AAᵧ::AminoAcid,
    substitutionmatrix::SubstitutionMatrix{AminoAcid,Int64}, similarityscorecutoff::Int64, similarityvalidator::Function
)
    similarityscore = substitutionmatrix[AAᵦ, AAᵧ]
    similarityvalidator(similarityscore, similarityscorecutoff)
end


"""
    anysimilarity(Sᵢ, Sⱼ, AA_groups)

Return `true` if at least one change from a `AAᵦ` to a `AAᵧ`, `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, is considered a similar change
according to `aagroups`' classification.   
That is, both `AAᵦ` and `AAᵧ` have the same classification and thus considered similar. 
"""
function anysimilarity(
    Sᵢ::Set{AminoAcid}, Sⱼ::Set{AminoAcid},
    aagroups::Dict{AminoAcid,String}
)
    ThreadsX.any(
        [
        issimilar(AAᵦ, AAᵧ, aagroups)
        for AAᵦ ∈ Sᵢ
        for AAᵧ ∈ Sⱼ
    ]
    )
end


"""
    anysimilarity(Sᵢ, Sⱼ, substitutionmatrix, similarityscorecutoff, similarityvalidator)
    
Return `true` if at least one change from a `AAᵦ` to a `AAᵧ`, `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, is considered a similar change
according to their substitution score in `substitutionmatrix`.   
That it, the substitution score is `>`/`≥`/`<`/`≤` `similarityscorecutoff`. 
The specific comparison (e.g., `≥`) is determined by `similarityvalidator`.
"""
function anysimilarity(
    Sᵢ::Set{AminoAcid}, Sⱼ::Set{AminoAcid},
    substitutionmatrix::SubstitutionMatrix{AminoAcid,Int64}, similarityscorecutoff::Int64, similarityvalidator::Function
)
    ThreadsX.any(
        [
        issimilar(AAᵦ, AAᵧ, substitutionmatrix, similarityscorecutoff, similarityvalidator)
        for AAᵦ ∈ Sᵢ
        for AAᵧ ∈ Sⱼ
    ]
    )
end