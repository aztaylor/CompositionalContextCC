using IJulia, DifferentialEquations, Plots

p = (
  Pᵣᵗ=18.931, #, [description="Total Pᵣ, in myTXTL no degredataion", connect=Flow, unit=u"nM"]
  pₗᵗ=11.0, #, [description="Total lac promoter", connect=Flow, unit=u"nM"]
  pₜᵗ=11.0, #, [description="Total Tet Promoter", connect=Flow, unit=u"nM"]
  lₗ=40, #, [description="Length of plac", unit=u"bp"]
  lₜ=44, #, [description="Length of ptet", unit=u"bp"]
  lₛ=240, #, [description="Length of mSpinach and T500 terminator", unit=u"bp"]
  lₘ=63, #,  [description="Length of MG and T500 terminator", unit=u"bp"]
  kᵢᵐ=7e-2, #, [description="Max initiation rate", unit=u"nM^-1*s^-1"]
  kₑᵐ=7e-2, #, [description="Max elongation rate", unit=u"s^-1"]
  kσₘₘ=50e-3, #, [description="MM constant for supercoiling hillfunctions",unit=u"nM"]
  kₒ=0.04, #, [description = "Rate of open complex formation", unit=u"s^-1"]
  kᵣ=1/170, #, [description = "RFP maturation rate", unit=u"s^-1"]
  kᵢₗₐ=6e3, #, [description="Rate of DNA-free apolacI IPTG binding", unit=u"nM^-1*s^-1"]
  kᵢₗᵤ=1, #,[description="Rate of apolacI IPTG disassociation", unit =u"s^-1"]
  kₚₗₐ=10, #, [description="lacI-promoter asossiation rate", unit=u"s^-1"]
  kₚₗᵤ=0.022, #, [description ="lacI-promoter disassociation rate", unit=u"s^-1"]
  kₐₜₐ=6e3, #, [description = "aTc-TetR association rate", unit=u"nM^-1*s^-1"]
  kₐₜᵤ=1, #, [description="aTc-TetR disassociation rate", unit=u"s^-1"]
  kₚₜₐ=10, #, [description="promoter-TetR association rate", unit=u"s^-1"]
  kₚₜᵤ =0.022, #, [description = "promoter-TetR disassociation rate", unit=u"s^-1"]
  ρₗ=0, #, [description="Rate of lacI production", unit=u"nM*s^-1"]
  ρₜ=0, #, [description="Rate of tetR production", unit=u"nM*s^-1"]
  δₛ=log(2)/(30*60), #, [description = "mSpinach degredation rate", unit=u"s^-1"]
  δₘ=log(2)/(60*60), #, [description = "MG degredation rate", unit=u"s^-1"]
  δₚ=0, #,[description="Average protein degredation rate", unit=u"s^-1"]
  σ₀=-0.065, #, [description="Natural B-form DNA supercoil state", unit=u"turn*bp^-1"]
  h₀=10.5 
)
    #Pᵣᵗ = Pᵣ+ecₛ+ecₘ+ccₗ+ccₜ
    #pₗᵗ = pₗ+ccₗ+ecₛ+pₗc
    #pₜᵗ = pₜ+ccₜ+ecₘ+pₜc
    #Iᵢᵗ = Iᵢ+aRₗ+pₗc
    #Iₐᵗ = Iₐ+aRₜ+pₜc
# Define auxilary functions that define and update the parameters wrt to the state of the 
# system
function n(σt₁, σt₋₁, p₋₁, p₋₁c, lₚ₁, lₜ₁, lₜ₋₁, h₀, lᵢ=150)
  """
  Computes the space between the promoter and the kink within the intergenic region.
  - σ: supercoiling state. P signifies the promoter region, t the ORF and terminator with 
    1 being the sense gene and -1 antisense. 0 indicates the natural  state of B-form DNA
  - h₀: The natural bp/turn of B-form DNA. 
  - p: Promoter conc. of the antisense-perspective gene. No c indicates free, 
    c indicates in complex with repressor.
  - l: Length of element. a (all) indicates linear DNA/plasmid, p promoter(1 is 
    sense-perspective, -1 antisense), t indicates ORF and terminator, i intergenic space.
  """
  Δₖᵢₙₖ = (σt₁+σt₋₁)*h₀
  n = lₚ₁+lₜ₁+(lᵢ/2)*p₋₁/(p₋₁+p₋₁c) + (lᵢ/2+lₜ₋₁)*p₋₁c/(p₋₁+p₋₁c)-Δₖᵢₙₖ
end 
function m(σ, σ₀=-0.065, T₀=2, G₀=18.93, τ=0.5, γ=0.5, kσₘₘ= 50e-3)
  """
  Accounts for changes in supercoiling due to Gyrase and Topoisomerase I.
  - σ: supercoiling state. 0 indicates the natural B-form DNA state.
  - T: Conc. of topoisomerase.
  - G: Conc. of gyrase.
  - τ: Rate topoisomerase introduces positive supercoils.
  - γ: Rate gyrase introduces negative supercoils.
  - kσₘₘ: Michaelis-menten constand associated with both enzymes.
  """
  σ⁺ = (σ+abs(σ))/2  
  σ⁻ = (σ-abs(σ))/2
  m=T₀*τ*(σ⁻)/kσₘₘ/(σ₀+abs(σ-σ₀))/kσₘₘ+G₀*γ*(σ⁺)/kσₘₘ/(σ₀+abs(σ-σ₀))/kσₘₘ
end

σᵖ(σ₀,l,lₚ=2829) = σ₀*(lₚ/l)
kᵢ(σ, σᵖ, kᵢₘ=7e-2) = σᵖ*kᵢₘ/(σᵖ+((σ-σᵖ)^2))
kₑ(σ, σᵖ, kₑₘ=7e-2) = σᵖ*kₑₘ/(σᵖ+((σ-σᵖ)^2))

# Define State Equations.
function reporterDynamics!(du, u, p, t)
# 1   2   3    4    5   6   7  8   9  10
  Cₛ, Cₘ, ecₛ, ecₘ, ccₗ, ccₜ, Rₗ, Rₜ, Iᵢ, Iₐ = u
  kₑₛ. kₑₘ, δₛ, δₘ, kₒ, kᵢₗ, kᵢₜ, kₑₛ, kₑₘ, kₚₗᵤ, kₚₜᵤ, kᵢₗᵤ, kₐₜᵤ, kᵢₗₐ, kₐₜₐ, kₑₘ, Pᵣᵗ = p
  du[1] = kₑₛ*ecₛ - δₛ*Cₛ
  du[2] = kₑₘ*ecₘ - δₘ*Cₘ
  du[3] = kₒ*ccₗ - kₑₛ*ecₛ
  du[4] = kₒ*ccₜ - kₑₘ*ecₘ
  du[5] = kᵢₗ*(Pᵣᵗ-ecₛ-ecₘ-ccₗ-ccₜ)*(pₗᵗ-ccₗ-ecₛ-pₗc) - (kₑₛ-kₒ)*ccₗ
  du[6] = kᵢₜ*(Pᵣᵗ-ecₘ-ecₛ-ccₗ-ccₜ)*(pₜᵗ-ccₜ-ecₘ-pₜc) - (kₑₘ-kₒ)*ccₜ
  du[7] = ρₗ + kᵢₗᵤ*(Iᵢᵗ-Iᵢ) + kₚₗᵤ*(Rₗᵗ-Rₗ-Iᵢᵗ-Iᵢ) - kᵢₗₐ*Rₗ*Iᵢ - kₚₗₐ*Rₗ - δₚ*Rₗ
  du[8] = ρₜ + kₐₜᵤ*(Iₐᵗ-Iₐ) + kₚₜᵤ*(Rₜᵗ-Rₜ-Iₐᵗ-Iₐ) - kₐₜₐ*Rₜ*Iₐ - kₚₜₐ*Rₜ - δₚ*Rₜ
  du[9] = kᵢₗₐ*(Rₗ+pₗc)*Iᵢ + kᵢₗᵤ*(Rₗᵗ-Rₗ-pₗc)
  du[10] = kₐₜₐ*(Rₜ+pₜc)*Iₐ + kₐₜᵤ*(Rₜᵗ-Rₜ-pₜc)
end
function convergentσDynaimcs!(du, u, p, t)
  σₜₛ, σₚₗ, σₜₘ, σₚₜ = u 
  lₗ, lₜ, lₛ, lₘ, δₛ, δₘ, h₀ = p
  du[11] = -(du[1]-δₛ*Cₛ-du[3])*(lₛ)/(2*h₀*nfₛ)-(du[3]-du[5])*(lₗ/2*h₀*nfₛ)+m
  du[12] = -(du[3]-du[5])*(lₗ/2*h₀*nfₛ)+m
  du[13] = -(du[2]-δₘ*Cₘ-du[4])*(lₘ)/(2*h₀*nfₘ)-(du[4]-du[6])*(lₜ/2*h₀*nfₘ)+m
  du[14] = -(du[4]-du[6])*(lₘ/2*h₀*nfₘ)+m
end 
function convergentSystem(du, u, p, t)
  reporterDynamics!(du[1:9], u[1:9], p, t)
  convergentσDynaimcs!(du[10:13], u[10:13], p, t)
end

sys =ODESystem(convergentSystem,t)
#@named sys = ODEProblem(sys,[], check_length=false)
