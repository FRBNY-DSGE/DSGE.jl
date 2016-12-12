#=
p = Dict(:z_star => "zstar", :rstar => "rstar", :Rstarn => "Rstarn", :r_k_star => "rkstar", :wstar => "wstar", :Lstar => "Lstar", :kstar => "kstar", :kbarstar => "kbarstar", :istar => "istar",
 :ystar => "ystar", :cstar => "cstar", :wl_c => "wl_c", :nstar => "nstar", :vstar => "vstar", :ζ_spσ_ω => "zeta_spsigw", :ζ_spμ_e => "zeta_spmue", :ζ_nRk => "zeta_nRk", :ζ_nR => "zeta_nR", :ζ_nqk => "zeta_nqk",
 :ζ_nn => "zeta_nn", :ζ_nμ_e => "zeta_nmue", :ζ_nσ_ω => "zeta_nsigw");
ss = Dict(p[v.key] => v.value for v in m.steady_state);
file = matopen("ss.mat", "w")
write(file, "ss", ss)
close(file)
=#
p=Dict(:α=>"alp",:ζ_p=>"zeta_p",:ι_p=>"iota_p",:δ=>"del",:Upsilon=>"ups",:Φ=>"Bigphi",:S′′=>"s2",:h=>"h",:ppsi=>"ppsi",:ν_l=>"nu_l",:ζ_w=>"zeta_w",:ι_w=>"iota_w",:λ_w=>"muwstar",:β=>"bet_",
:ψ1=>"psi1",:ψ2=>"psi2",:ψ3=>"psi3",:π_star=>"pistar_",:σ_c=>"sigmac",:ρ=>"rho",:ϵ_p=>"epsp",:ϵ_w=>"epsw",:Fω=>"Fom_",:spr=>"sprd_",:ζ_spb=>"zeta_spb",:γ_star=>"gammstar",:γ=>"gam_",:Lmean=>"Lmean",
:g_star=>"gstar",:ρ_g=>"rho_g",:ρ_b=>"rho_b",:ρ_μ=>"rho_mu",:ρ_z=>"rho_z",:ρ_λ_f=>"rho_laf",:ρ_λ_w=>"rho_law",:ρ_rm=>"rho_rm",:ρ_σ_w=>"rho_sigw",:ρ_μ_e=>"rho_mue",:ρ_γ=>"rho_gamm",:ρ_π_star=>"rho_pist",:ρ_lr=>"rho_lr",:ρ_z_p=>"rho_zp",
:ρ_tfp=>"rho_tfp",:ρ_gdpdef=>"rho_gdpdef",:ρ_corepce=>"rho_pce",:σ_g=>"std_g_sh",:σ_b=>"std_b_sh",:σ_μ=>"std_mu_sh",:σ_z=>"std_z_sh",:σ_λ_f=>"std_laf_sh",:σ_λ_w=>"std_law_sh",:σ_r_m=>"std_rm_sh",:σ_σ_ω=>"std_sigw_sh",:σ_μ_e=>"std_mue_sh",:σ_γ=>"std_gamm_sh",:σ_π_star=>"std_pist_sh",
:σ_lr=>"std_lr_sh",:σ_z_p=>"std_zp_sh",:σ_tfp=>"std_tfp_sh",:σ_gdpdef=>"std_gdpdef_sh",:σ_corepce=>"std_pce_sh",
:η_gz=>"eta_gz",:η_λ_f=>"eta_laf",:η_λ_w=>"eta_law",:Iendoα=>"modelalp_ind",:Γ_gdpdef=>"gamm_gdpdef",:δ_gdpdef=>"del_gdpdef")
p=merge(p,Dict(Symbol("σ_r_m$i")=>"std_rm_sh$i" for i = 1:20));
#=
params = Dict(p[θ.key] => θ.value for θ in m.parameters);
file = matopen("params.mat", "w")
write(file, "params", params)
close(file)
=#
pr=Dict()
for θ in m.parameters
    if !θ.fixed
      pr[p[θ.key]] = logpdf(θ)
    end
end
println(pr)
file = matopen("prior.mat", "w")
write(file, "prior", pr)
close(file)
readline(STDIN)
quit()
