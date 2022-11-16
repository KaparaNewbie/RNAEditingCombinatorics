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