using Catalyst
Retina = @reaction_network begin
    µ, 0 --> Sci0
    φ, 0 --> Sci1
    β, 0 --> Sci2 
    θ, Ce0 --> Ce0 + C0 
    λ, Ce1 --> Ce1 + C1
    κ, Ce2 --> Ce2 + C2 
    η, Sci0 + C0 --> 0
    ψ, Sci1 + C1 --> 0
    ɛ̝, Sci2 + C2 --> 0
    k1, Sci0 --> Sci0 + Sr0 + Sh0 + Sp0 + Sb0 
    p1, Sci1 --> Sci1 + Sb1 + Sh1 
    p26, Sci2 --> Sci2 + Sh2 + Sr2 + Sp2
    k17, Sr0 --> Sr0 + R0 
    p33, Sr2 --> Sr2 + R2
    k18, Sh0 --> Sh0 + HR0 
    p7, Sh1 --> Sh1 + HR1
    p19, Sh2 --> Sh2 + HR2
    k19, Sp0 --> Sp0 + P0 
    p37, Sp2 --> Sp2 + P2
    k28, Sb0 --> Sb0 + Rb0  
    p9, Sb1 --> Sb1 + Rb1 
    k2, Sr0 --> 0
    p32, Sr2 --> 0
    k10, Sh0 --> 0
    p3, Sh1 --> 0
    p18, Sh2 --> 0
    k13, Sp0 --> 0
    p31, Sp2 --> 0
    k27, Sb0 --> 0
    p2, Sb1 --> 0
    k3, R0 + CL0 --> R0 + CL0 + Cf0 
    p34, R2 + CL2 --> R2 + Cf2 
    k26, Rb0 + CH0 --> Rb0 + CH0 + Cf0 
    p25, Rb1 + CH1 --> Rb1 + Cf1 
    k11, R0 --> 0
    p35, R2 --> 0
    k29, Rb0 --> 0
    p10, Rb1 --> 0
    k14, P0 + R0 --> P0 
    k16, P0 --> 0
    p39, P2 + R2 --> P2 
    p38, P2 --> 0
    k4, Cf0 --> Cp0 
    p4, Cf1 --> Cp1 
    p36, Cf2 --> Cp2 
    p11, Cp1 --> CO 
    (k5, k6), Cp0 <--> Ce0 
    (p5, p6), Cp1 <--> Ce1 
    (p28, p29), Cp2 <--> Ce2 
    (k7, k8), Ce0 <--> E0 
    (p16, p17), Ce1 <--> E1 
    α, 0 --> H0 
    ω, 0 --> H1 
    δ, 0 --> H2 
    k9, H0 + HR0 --> HR0 + Ce0 
    p12, H1 + HR1 --> HR1 + Ce1 
    p21, H2 + HR2 --> HR2 + Ce2 
    k15, Ce0 + HR0 --> Ce0 
    p20, Ce2 + HR2 --> Ce2
    k12, Cp0 + CO --> Cg 
    k21, Ce0 --> 0
    p15, Ce1 --> 0
    p13, Ce2 --> 0
    k22, Ce0 + Cg --> B 
    k24, B --> 0
end
