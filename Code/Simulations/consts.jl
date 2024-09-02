const ∅ = Set()

# https://ib.bioninja.com.au/standard-level/topic-2-molecular-biology/24-proteins/amino-acids.html 
const AA_groups = Dict(
    AA_G => "NON-POLAR",  # Gly
    AA_A => "NON-POLAR",  # Ala
    AA_V => "NON-POLAR",  # Val
    AA_C => "NON-POLAR",  # Cys
    AA_P => "NON-POLAR",  # Pro
    AA_L => "NON-POLAR",  # Leu
    AA_I => "NON-POLAR",  # Ile
    AA_M => "NON-POLAR",  # Met
    AA_W => "NON-POLAR",  # Trp
    AA_F => "NON-POLAR",  # Phe

    AA_S => "POLAR",  # Ser
    AA_T => "POLAR",  # Thr
    AA_Y => "POLAR",  # Tyr
    AA_N => "POLAR",  # Asn
    AA_Q => "POLAR",  # Gln

    AA_K => "POSITIVE",  # Lys
    AA_R => "POSITIVE",  # Arg
    AA_H => "POSITIVE",  # His

    AA_D => "NEGATIVE",  # Asp
    AA_E => "NEGATIVE",  # Glu
)




"""
Classification by volume and polarity (Miyata et al. 1979).  

https://www.nature.com/articles/nmicrobiol201791#Sec7 - Supplementary Table 1  
"""
const AA_groups_Miyata1979 = Dict(

    AA_C => "Special",  # Cys
    
    AA_A => "Neutral-Small",  # Ala
    AA_G => "Neutral-Small",  # Gly
    AA_P => "Neutral-Small",  # Pro
    AA_S => "Neutral-Small",  # Ser
    AA_T => "Neutral-Small",  # Thr
    
    AA_N => "Polar-RelativelySmall",  # Asn
    AA_D => "Polar-RelativelySmall",  # Asp
    AA_Q => "Polar-RelativelySmall",  # Gln
    AA_E => "Polar-RelativelySmall",  # Glu

    AA_R => "Polar-RelativelyLarge",  # Arg
    AA_H => "Polar-RelativelyLarge",  # His
    AA_K => "Polar-RelativelyLarge",  # Lys

    AA_I => "NonPolar-RelativelySmall",  # Ile
    AA_L => "NonPolar-RelativelySmall",  # Leu
    AA_M => "NonPolar-RelativelySmall",  # Met
    AA_V => "NonPolar-RelativelySmall",  # Val
    
    AA_F => "NonPolar-RelativelyLarge",  # Phe
    AA_W => "NonPolar-RelativelyLarge",  # Trp
    AA_Y => "NonPolar-RelativelyLarge",  # Tyr 

)