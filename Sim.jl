using IJulia, ModelingToolkit, DifferentialEquations, Plots

@variables t #, [unit=u"s"]
D = Differential(t)

ConservationLaws = @connector ConsLaws begin
  @parameters begin
    Pᵣᵗ=18.931e3 #, [description="Total Pᵣ, in myTXTL no degredataion", connect=Flow, unit=u"nM"]
    pₗᵗ=11.0 #, [description="Total lac promoter", connect=Flow, unit=u"nM"]
    pₜᵗ=11.0 #, [description="Total Tet Promoter", connect=Flow, unit=u"nM"]
    Rₗᵗ=10 #, [description="Total Lac Repressor, in myTXTL no degredataion", connect=Flow, unit=u"nM"]
    Rₜᵗ=0.0 #, [description="Total Tet Repressor, in myTXTL no degredataion", connect=Flow, unit=u"nM"]
    Iₐᵗ=0.0 #, [description="Total IPTG, not metabolized by the reaction volume", connect=Flow, unit=u"nM"]
    Iᵢᵗ=0.0 #, [description="Total aTc, not metabolized by the reaction volume", connect=Flow, unit=u"nM"]
  end
  @variables (begin
    ecₛ(t)=0, [input],#, [description="Number of mSpinach elongation comlexes", unit=u"nM"]
    ecₘ(t)=0, [input], #, [description="Number of mSpinach closed dna comlexes", unit=u"nM"]
    ccₛ(t)=0, [input], #, [description="Number of MG elongation comlexes", unit=u"nM"]
    ccₘ(t)=0, [input],#, [description="Number of MG closed dna comlexes", unit=u"nM"]
    pₗ(t)=pₗᵗ, [connect = Flow], #, [description="conc plac", unit=u"nM"]
    pₜ(t)=pₜᵗ, [connect = Flow], #, [description="conc pTet", unit=u"nM"]
    pₗc(t)=0, [input],  #, [description="conc plac-lacI complex", unit=u"nM"]
    pₜc(t)=0, [input], #, [description="conc pTet-TetR complex", unit=u"nM"]
    Rₗ(t)=0, [input], #, [description="conc LacI Repressor", unit=u"nM"]
    Rₜ(t)=0, [input], #, [description="conc TetR Repressor", unit=u"nM"]   
    aRₗ(t)=0, [input], #, [description="conc apo LacI", unit=u"nM"]
    aRₜ(t)=0, [input], #, [description="conc apo TetR", unit=u"nM"]
    Iᵢ(t)=Iᵢᵗ, [connect = Flow], #, [description="conc IPTG", unit=u"nM"]
    Iₐ(t)=Iₐᵗ, [connect = Flow],#, [description="conc aTc", unit=u"nM"]
    Pᵣ(t)=Pᵣᵗ, [connect = Flow] #, [description="conc Pᵣ", unit=u"nM"]
  end)
  @equations begin
    Pᵣᵗ ~ Pᵣ+ecₛ+ecₘ+ccₛ+ccₘ
    pₗᵗ ~ pₗ+ccₛ+ecₛ+pₗc
    pₜᵗ ~ pₜ+ccₘ+ecₘ+pₜc
    Iᵢᵗ ~ Iᵢ+aRₗ+pₗc
    Iₐᵗ ~ Iₐ+aRₜ+pₜc
  end
end

RateLaws = @connector rates begin
  @parameters begin
    lₗ=40 #, [description="Length of plac", unit=u"bp"]
    lₜ=44 #, [description="Length of ptet", unit=u"bp"]
    lₛ=240 #, [description="Length of mSpinach and T500 terminator", unit=u"bp"]
    lₘ=63 #,  [description="Length of MG and T500 terminator", unit=u"bp"]
    lₚ=2892 #, [description="Length of plasmid", unit=u"bp"]
    kᵢₘₐₓ=7e-2 #, [description="Max initiation rate", unit=u"nM^-1*s^-1"]
    kₑₘₐₓ=7e-2 #, [description="Max elongation rate", unit=u"s^-1"]
    kσₘₘ=50e-3 #, [description="MM constant for supercoiling hillfunctions",unit=u"nM"]
    kₒ=0.04 #, [description = "Rate of open complex formation", unit=u"s^-1"]
    kᵣ=1/170 #, [description = "RFP maturation rate", unit=u"s^-1"]
    kᵢₗₐ=6e3 #, [description="Rate of DNA-free apolacI IPTG binding", unit=u"nM^-1*s^-1"]
    kᵢₗᵤ=1 #,[description="Rate of apolacI IPTG disassociation", unit =u"s^-1"]
    kₚₗₐ=10 #, [description="lacI-promoter asossiation rate", unit=u"s^-1"]
    kₚₗᵤ=0.022 #, [description ="lacI-promoter disassociation rate", unit=u"s^-1"]
    kₐₜₐ=6e3 #, [description = "aTc-TetR association rate", unit=u"nM^-1*s^-1"]
    kₐₜᵤ=1 #, [description="aTc-TetR disassociation rate", unit=u"s^-1"]
    kₚₜₐ=10 #, [description="promoter-TetR association rate", unit=u"s^-1"]
    kₚₜᵤ =0.022 #, [description = "promoter-TetR disassociation rate", unit=u"s^-1"]
    ρₗ=0 #, [description="Rate of lacI production", unit=u"nM*s^-1"]
    ρₜ=0 #, [description="Rate of tetR production", unit=u"nM*s^-1"]
    δₛ=log(2)/(30*60) #, [description = "mSpinach degredation rate", unit=u"s^-1"]
    δₘ=log(2)/(60*60) #, [description = "MG degredation rate", unit=u"s^-1"]
    δₚ=0 #,[description="Average protein degredation rate", unit=u"s^-1"]
    σ₀=-0.065 #, [description="Natural B-form DNA supercoil state", unit=u"turn*bp^-1"] 
    σᵒpₗ=σ₀*lₚ/lₗ #, [description="Approximate Optimal supercoiling density, plac"]
    σᵒtₛ=σ₀*lₚ/lₛ #, [description="Approximate Optimal supercoiling density, plac"]
    σᵒpₜ=σ₀*lₚ/lₜ #, [description="Approximate Optimal supercoiling density, pTet"]
    σᵒtₘ=σ₀*lₚ/lₘ #, [description="Approximate Optimal supercoiling density, pTet"]
  end
  @variables begin
    σtₛ(t)=-6, [input],#, [description="supercoil state of mSpinach ORF"]
    σtₘ(t)=-3, [input], #, [description="supercoil state of MG ORF"]
    σpₗ(t)=-6, [input], #, [description="supercoil state of mSpinach promoter"]
    σpₜ(t)=-3, [input], #, [description="supercoil state of MG promoter"]
    kᵢₗ(t)=0.5, [output], #, [description="transcription initiation rate plac", unit=u"nM^-1*s^-1"]
    kᵢₜ(t)=0.5, [output], #, [description="transcription initiation rate pTet", unit=u"nM^-1*s^-1"]
    kₑₛ(t)=0.5,[output], #, [description="transcription elongation rate mSpinach", unit=u"s^-1"]
    kₑₘ(t)=0.5, [output] #, [description="transcription elongation rate MG", unit=u"s^-1"]
  end
  @equations begin
    kᵢₗ ~ σᵒpₗ*kᵢₘₐₓ/(σᵒpₗ+((σpₗ-σᵒpₗ)^2))
    kₑₛ ~ σᵒtₛ*kₑₘₐₓ/(σᵒtₛ +((σtₛ-σᵒtₛ)^2))
    kᵢₜ ~ σᵒpₜ*kᵢₘₐₓ/(σᵒpₜ+((σpₜ-σᵒpₜ)^2))
    kₑₘ ~ σᵒtₘ*kₑₘₐₓ/(σᵒtₘ+((σtₘ-σᵒtₘ)^2))
  end
end

ReporterLaws = @mtkmodel reporterDynamics begin
  @components begin
    c = ConsLaws()
    r = rates()
  end

  @variables begin
    Cₛ(t)=0 #, [description="mSpinach Transcript", unit=u"nM"]
    Cₘ(t)=0 #, [description="MG Transcript", unit=u"nM"]
  end
  @equations begin
    D(Cₛ) ~ r.kₑₛ*c.ecₛ - r.δₛ*c.Cₛ
    D(Cₘ) ~ r.kₑₘ*c.ecₘ - r.δₘ*c.Cₘ
    D(ecₛ) ~ r.kₒ*c.ccₛ - r.kₑₛ*c.ecₛ
    D(ecₘ) ~ r.kₒ*c.ccₘ - r.kₑₘ*c.ecₘ
    D(ccₛ) ~ r.kᵢₗ*(c.Pᵣᵗ-c.ecₛ-c.ecₘ-c.ccₛ-c.ccₘ)*(c.pₗᵗ-c.ccₛ-c.ecₛ-c.pₗc) - (r.kₑₛ-r.kₒ)*c.ccₛ
    D(ccₘ) ~ r.kᵢₜ*(c.Pᵣᵗ-c.ecₘ-c.ecₛ-c.ccₛ-c.ccₘ)*(c.pₜᵗ-c.ccₘ-c.ecₘ-c.pₜc) - (r.kₑₘ-r.kₒ)*c.ccₘ
    D(Rₗ) ~ r.ρₗ + r.kᵢₗᵤ*(c.Iᵢᵗ-c.Iᵢ) + r.kₚₗᵤ*(c.Rₗᵗ-c.Rₗ-c.Iᵢᵗ-c.Iᵢ) - r.kᵢₗₐ*c.Rₗ*c.Iᵢ - r.kₚₗₐ*c.Rₗ - r.δₚ*c.Rₗ
    D(Rₜ) ~ r.ρₜ + r.kₐₜᵤ*(c.Iₐᵗ-c.Iₐ) + r.kₚₜᵤ*(c.Rₜᵗ-c.Rₜ-c.Iₐᵗ-c.Iₐ) - r.kₐₜₐ*c.Rₜ*c.Iₐ - r.kₚₜₐ*c.Rₜ - r.δₚ*c.Rₜ
    D(Iᵢ) ~ r.kᵢₗₐ*(c.Rₗ+c.pₗc)*c.Iᵢ + r.kᵢₗᵤ*(c.Rₗᵗ-c.Rₗ-c.pₗc)
    D(Iₐ) ~ r.kₐₜₐ*(c.Rₜ+c.pₜc)*c.Iₐ + r.kₐₜᵤ*(c.Rₜᵗ-c.Rₜ-c.pₜc)
  end
end

σMaintenanceLaws = @mtkmodel σMaintenanceDynamics begin
  @extend σDynamics()
  @parameters begin
    G₀=18.93 #, [description="Concentration Gyrase", unit=u"μM"]
    T₀=2 #, [description="Conc Topoisomerase", unit=u"μM"]
    τ=0.5 #, [description="Rate of topoisomerase activity", unit=u"turn*s^-1"]
    γ=0.5 #, [description="Rate of Gyrase activit", unit=u"turn*s^-1"]
    kσₘₘ=200 #, [description="Michaelis-Menten constant for gyrase", unit=u"μM"]
  end 
  @variables begin
    σpₗ⁺(t) = 0 #, [description="Decomposition of the plac supercoiling state into strictly positive parts"] 
    σpₗ⁻(t) = 0 #, [description="Decomposition of the plac supercoiling state into strictly negative parts"]
    σtₛ⁺(t) = 0 #, [description="Decomposition of the mSpinach supercoiling state into strictly positive parts"]
    σtₛ⁻(t) = 0 #, [description="Decomposition of the mSpinach supercoiling state into strictly negative parts"]
    σpₜ⁺(t) = 0 #, [description="Decomposition of the pTet supercoiling state into strictly positive parts"]
    σpₜ⁻(t) = 0 #, [description="Decomposition of the pTet supercoiling state into strictly negative parts"]
    σtₘ⁺(t) = 0 #, [description="Decomposition of the MG supercoiling state into strictly positive parts"]
    σtₘ⁻(t) = 0 #, [description="Decomposition of the MG supercoiling state into strictly negative parts"]
    mpₗ(t) = 0 #, [description="Maintenance Dynamics from topoisomerase and Gyrase", unit=u"nM*s^-1"]
    mtₛ(t) = 0 #, [description="Maintenance Dynamics from topoisomerase and Gyrase", unit=u"nM*s^-1"]
    mpₜ(t) = 0 #, [description="Maintenance Dynamics from topoisomerase and Gyrase", unit=u"nM*s^-1"]
    mtₘ(t) = 0 #, [description="Maintenance Dynamics from topoisomerase and Gyrase", unit=u"nM*s^-1"]
  end  
  @equations begin
    σpₗ⁺ = (σpₗ+abs(σpₗ))/2  
    σpₗ⁻ = (σpₗ-abs(σpₗ))/2
    σtₛ⁺ = (σtₛ+abs(σtₛ))/2
    σtₛ⁻ = (σtₛ-abs(σtₛ))/2
    σpₜ⁺ = (σpₜ+abs(σpₜ))/2
    σpₜ⁻ = (σpₜ-abs(σpₜ))/2
    σtₘ⁺ = (σtₘ+abs(σtₘ))/2
    σtₘ⁻ = (σtₘ-abs(σtₘ))/2
    mpₗ~T₀*τ*(σpₗ⁻)/kσₘₘ/(σ₀+abs(σpₗ-σ₀))/kσₘₘ+G₀*γ*(σpₗ⁺)/kσₘₘ/(σ₀+abs(σpₗ-σ₀))/kσₘₘ
    mtₛ~T₀*τ*(σtₛ⁻)/kσₘₘ/(σ₀+abs(σtₛ-σ₀))/kσₘₘ+G₀*γ*(σtₛ⁺)/kσₘₘ/(σ₀+abs(σtₛ-σ₀))/kσₘₘ
    mpₜ~T₀*τ*(σpₜ⁻)/kσₘₘ/(σ₀+abs(σpₜ-σ₀))/kσₘₘ+G₀*γ*(σpₜ⁺)/kσₘₘ/(σ₀+abs(σpₜ-σ₀))/kσₘₘ
    mtₘ~T₀*τ*(σtₘ⁻)/kσₘₘ/(σ₀+abs(σtₘ-σ₀))\kσₘₘ+G₀*γ*(σtₘ⁺)/kσₘₘ/abs(σ₀+(σtₘ-σ₀))/kσₘₘ
  end
end

interspaceLaws = @mtkmodel nᵢ begin
  @extend σDyanmicsMaintanence
  @parameters begin
    h₀= 10.5 #, [description="basepairs per right-hand turn", unit=u"bp*turn^-1"]
    σ₀=-0.065 #, [description="standard supercoil state"]
    lᵢ = 150 #, [description="Intergenic Spacer Length", unit="bp"]
  end
  @variables begin  
    nfₛ(t)=100 #, [description="length of promoter, mSpinach ORF, terminator and available intergenic space", unit=u"bp"]
    nfₘ(t)=100 #, [description="length of promoter, MG ORF, terminator and available intergenic space", unit=u"bp"]
  end
  @equations begin
    Δₖᵢₙₖ = (r.σtₛ+r.σtₘ)*h₀
    nfₛ = (c.lₗ+c.lₛ+lᵢ/2(c.pₜ/(c.pₜ+c.pₜc))+(lᵢ/2+c.lₘ)*c.pₜc/(c.pₜ+c.pₜc)-Δₖᵢₙₖ) 
    nfₘ =(c.lₜ+c.lₘ+lᵢ/2(c.pₗ/(c.pₗ-c.pₗc))+(lᵢ/2+c.lₛ)*c.pₗc/(c.pₗ+c.pₗc)-Δₖᵢₙₖ)
  end
end

σLaws = @mtkmodel σDynamics begin
  @extend nᵢ
  @equations begin
    D(σtₛ) ~ -(D(Cₛ)-δₛ*Cₛ-D(ecₛ))*(lₛ)/(2*h₀*nfₛ)-(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+mtₛ
    D(σpₗ) ~ -(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+mpₗ
    D(σtₘ) ~ -(Cₘ-δₛ*Cₘ-D(ecₘ))*(lₘ)/(2*h₀*nfₘ)-(D(ecₘ)-D(ccₘ))*(lₜ/2*h₀*nfₘ)+mpₘ
    D(σpₜ) ~ -(D(ecₘ)-D(ccₘ))*(lₘ/2*h₀*nₘ)+mpₜ
  end
end


