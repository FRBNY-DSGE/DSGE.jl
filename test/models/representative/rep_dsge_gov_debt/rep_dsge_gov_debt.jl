using DSGE
using Test, BenchmarkTools
using JLD2

# What do you want to do?
check_steady_state = true
check_jacobian     = true
check_solution     = true
check_irfs         = true

path = dirname(@__FILE__)

het = HetDSGEGovDebt(testing_gamma = true, ref_dir = HETDSGEGOVDEBT)
het <= Setting(:steady_state_only, true)
steadystate!(het)

m   = RepDSGEGovDebt(het)

# Steady-state computation
if check_steady_state
    steadystate!(m, het)

    file = jldopen("$path/../../heterogeneous_agent/het_dsge_gov_debt/reference/steady_state.jld2", "r")
    saved_Rk   = read(file, "Rk")
    saved_ω    = read(file, "omega")
    saved_kl   = read(file, "kl")
    saved_k    = read(file, "k")
    saved_y    = read(file, "y")
    saved_T    = read(file, "T")

    close(file)

    file = jldopen("$path/reference/steady_state.jld2", "r")

    saved_c    = read(file, "c")
    saved_ell  = read(file, "l")
    saved_β    = read(file, "beta")

    close(file)

    @testset "Check steady state outputs" begin
        @test saved_Rk ≈ m[:Rkstar].value
        @test saved_ω  ≈ m[:ωstar].value
        @test saved_kl ≈ m[:klstar].value
        @test saved_k  ≈ m[:kstar].value
        @test saved_y  ≈ m[:ystar].value
        @test saved_T  ≈ m[:Tstar].value

        @test saved_c   ≈ m[:cstar].value
        @test saved_ell ≈ m[:lstar].value
        @test saved_β   ≈ m[:βstar].value
    end
end

# Jacobian computation
if check_jacobian
    steadystate!(m, het)
    JJ = DSGE.jacobian(m)

    file = jldopen("$path/reference/jacobian.jld2", "r")
    saved_JJ = read(file, "JJ")

    nvars = read(file, "nvars")

    KKP  = read(file, "KKP")
    LRRP = read(file, "LRRP")
    LIIP = read(file, "LIIP")
    LYP = read(file, "LYP")
    LWP = read(file, "LWP")
    LXP = read(file, "LXP")
    BGP = read(file, "BGP")
    BP = read(file, "BP")
    GP = read(file, "GP")
    ZP = read(file, "ZP")
    MUP = read(file, "MUP")
    LAMWP = read(file, "LAMWP")
    LAMFP = read(file, "LAMFP")
    MONP = read(file, "MONP")
    RRP = read(file, "RRP")
    IIP = read(file, "IIP")
    TTP = read(file, "TTP")
    WP = read(file, "WP")
    HHP = read(file, "HHP")
    PIP = read(file, "PIP")
    PIWP = read(file, "PIWP")
    LAMP = read(file, "LAMP")
    YP = read(file, "YP")
    XP = read(file, "XP")
    MCP = read(file, "MCP")
    QP = read(file, "QP")
    RKP = read(file, "RKP")
    ELLP = read(file, "ELLP")
    TGP = read(file, "TGP")

    KK = KKP + nvars
    LRR = LRRP + nvars
    LII = LIIP + nvars
    LY = LYP + nvars
    LW = LWP + nvars
    LX = LXP + nvars
    BG = BGP + nvars
    B = BP + nvars
    G = GP + nvars
    Z = ZP + nvars
    MU = MUP + nvars
    LAMW = LAMWP + nvars
    LAMF = LAMFP + nvars
    MON = MONP + nvars
    RR = RRP + nvars
    II = IIP + nvars
    TT = TTP + nvars
    W = WP + nvars
    HH = HHP + nvars
    PI = PIP + nvars
    PIW = PIWP + nvars
    LAM = LAMP + nvars
    Y = YP + nvars
    X = XP + nvars
    MC = MCP + nvars
    Q = QP + nvars
    RK = RKP + nvars
    ELL = ELLP + nvars
    TG = TGP + nvars

    F1 = read(file, "F1")
    F5 = read(file, "F5")
    F6 = read(file, "F6")
    F7 = read(file, "F7")
    F8 = read(file, "F8")
    F9 = read(file, "F9")
    F10 = read(file, "F10")
    F11 = read(file, "F11")
    F12 = read(file, "F12")
    F13 = read(file, "F13")
    F14 = read(file, "F14")
    F15 = read(file, "F15")
    F16 = read(file, "F16")
    F17 = read(file, "F17")
    F18 = read(file, "F18")
    F40 = read(file, "F40")
    F41 = read(file, "F41")
    F19 = read(file, "F19")
    F20 = read(file, "F20")
    F24 = read(file, "F24")
    F31 = read(file, "F31")
    F32 = read(file, "F32")
    F33 = read(file, "F33")
    F34 = read(file, "F34")
    F35 = read(file, "F35")
    F36 = read(file, "F36")
    F37 = read(file, "F37")
    F38 = read(file, "F38")
    F39 = read(file, "F39")
    close(file)

    @testset "Euler equation" begin
        @test saved_JJ[F1,ELLP] ≈ JJ[F1, ELLP]
        @test saved_JJ[F1,ZP]   ≈ JJ[F1, ZP]
        @test saved_JJ[F1,B]    ≈ JJ[F1, B]
        @test saved_JJ[F1,ELL]  ≈ JJ[F1, ELL]
        @test saved_JJ[F1,RR]   ≈ JJ[F1, RR]
    end

    @testset "mkt clearing" begin
        @test saved_JJ[F5,Y]   ≈ JJ[F5, Y]
        @test saved_JJ[F5,G]   ≈ JJ[F5, G]
        @test saved_JJ[F5,X]   ≈ JJ[F5, X]
        @test saved_JJ[F5,ELL] ≈ JJ[F5, ELL]
    end

    @testset "lambda = average marginal utility" begin
        @test saved_JJ[F6,LAM] ≈ JJ[F6, LAM]
        @test saved_JJ[F6,ELL] ≈ JJ[F6, ELL]
    end

    @testset "transfer" begin
        @test saved_JJ[F7,TT] ≈ JJ[F7,TT]
        @test saved_JJ[F7,RK] ≈ JJ[F7,RK]
        @test saved_JJ[F7,KK] ≈ JJ[F7,KK]
        @test saved_JJ[F7,Z]  ≈ JJ[F7,Z]
        @test saved_JJ[F7,X]  ≈ JJ[F7,X]
        @test saved_JJ[F7,MC] ≈ JJ[F7,MC]
        @test saved_JJ[F7,TG] ≈ JJ[F7,TG]
    end

    @testset "investment" begin
        @test saved_JJ[F8,Q]  ≈ JJ[F8,Q]
        @test saved_JJ[F8,MU] ≈ JJ[F8,MU]
        @test saved_JJ[F8,XP] ≈ JJ[F8,XP]
        @test saved_JJ[F8,ZP] ≈ JJ[F8,ZP]
        @test saved_JJ[F8,X]  ≈ JJ[F8,X]
        @test saved_JJ[F8,LX] ≈ JJ[F8,LX]
        @test saved_JJ[F8,Z]  ≈ JJ[F8,Z]
    end

    @testset "tobin's q" begin
        @test saved_JJ[F9,LAM]  ≈ JJ[F9,LAM]
        @test saved_JJ[F9,LAMP] ≈ JJ[F9,LAMP]
        @test saved_JJ[F9,Q]    ≈ JJ[F9,Q]
        @test saved_JJ[F9,ZP]   ≈ JJ[F9,ZP]
        @test saved_JJ[F9,RKP]  ≈ JJ[F9,RKP]
        @test saved_JJ[F9,QP]   ≈ JJ[F9,QP]
    end

    @testset "capital accumulation" begin
        @test saved_JJ[F10,KKP] ≈ JJ[F10,KKP]
        @test saved_JJ[F10,KK]  ≈ JJ[F10,KK]
        @test saved_JJ[F10,Z]   ≈ JJ[F10,Z]
        @test saved_JJ[F10,MU]  ≈ JJ[F10,MU]
        @test saved_JJ[F10,X]   ≈ JJ[F10,X]
    end

    @testset "wage phillips curve" begin
        @test saved_JJ[F11,PIW]  ≈ JJ[F11,PIW]
        @test saved_JJ[F11,LAMW] ≈ JJ[F11,LAMW]
        @test saved_JJ[F11,HH]   ≈ JJ[F11,HH]
        @test saved_JJ[F11,W]    ≈ JJ[F11,W]
        @test saved_JJ[F11,PIWP] ≈ JJ[F11,PIWP]
        @test saved_JJ[F11,LAM]  ≈ JJ[F11,LAM]
    end

    @testset "price phillips curve" begin
        @test saved_JJ[F12,PI]   ≈ JJ[F12,PI]
        @test saved_JJ[F12,MC]   ≈ JJ[F12,MC]
        @test saved_JJ[F12,LAMF] ≈ JJ[F12,LAMF]
        @test saved_JJ[F12,PIP]  ≈ JJ[F12,PIP]
    end

    @testset "marginal cost" begin
        @test saved_JJ[F13,MC] ≈ JJ[F13,MC]
        @test saved_JJ[F13,W]  ≈ JJ[F13,W]
        @test saved_JJ[F13,RK] ≈ JJ[F13,RK]
    end

    @testset "gdp" begin
        @test saved_JJ[F14,Y]  ≈ JJ[F14,Y]
        @test saved_JJ[F14,Z]  ≈ JJ[F14,Z]
        @test saved_JJ[F14,KK] ≈ JJ[F14,KK]
        @test saved_JJ[F14,HH] ≈ JJ[F14,HH]
    end

    @testset "optimal k/l ratio" begin
        @test saved_JJ[F15,RK] ≈ JJ[F15,RK]
        @test saved_JJ[F15,W]  ≈ JJ[F15,W]
        @test saved_JJ[F15,HH] ≈ JJ[F15,HH]
        @test saved_JJ[F15,KK] ≈ JJ[F15,KK]
        @test saved_JJ[F15,Z]  ≈ JJ[F15,Z]
    end

    @testset "taylor rule" begin
        @test saved_JJ[F16,II]  ≈ JJ[F16,II]
        @test saved_JJ[F16,LII] ≈ JJ[F16,LII]
        @test saved_JJ[F16,PI]  ≈ JJ[F16,PI]
        @test saved_JJ[F16,Y]   ≈ JJ[F16,Y]
        @test saved_JJ[F16,LY]  ≈ JJ[F16,LY]
        @test saved_JJ[F16,Z]   ≈ JJ[F16,Z]
        @test saved_JJ[F16,MON] ≈ JJ[F16,MON]
    end

    @testset "fisher eqn" begin
        @test saved_JJ[F17,RR]  ≈ JJ[F17,RR]
        @test saved_JJ[F17,PIP] ≈ JJ[F17,PIP]
        @test saved_JJ[F17,II]  ≈ JJ[F17,II]
    end

    @testset "wage inflation" begin
        @test saved_JJ[F18,PIW] ≈ JJ[F18,PIW]
        @test saved_JJ[F18,PI]  ≈ JJ[F18,PI]
        @test saved_JJ[F18,Z]   ≈ JJ[F18,Z]
        @test saved_JJ[F18,W]   ≈ JJ[F18,W]
        @test saved_JJ[F18,LW]  ≈ JJ[F18,LW]
    end

    @testset "fiscal rule" begin
        @test saved_JJ[F40,TG]  ≈ JJ[F40,TG]
        @test saved_JJ[F40,RR]  ≈ JJ[F40,RR]
        @test saved_JJ[F40,BG]  ≈ JJ[F40,BG]
        @test saved_JJ[F40,Z]   ≈ JJ[F40,Z]
        @test saved_JJ[F40,Y]   ≈ JJ[F40,Y]
        @test saved_JJ[F40,G]   ≈ JJ[F40,G]
    end

    @testset "govt budget constraint" begin
        @test saved_JJ[F41,BGP] ≈ JJ[F41,BGP]
        @test saved_JJ[F41,RR]  ≈ JJ[F41,RR]
        @test saved_JJ[F41,BG]  ≈ JJ[F41,BG]
        @test saved_JJ[F41,Z]   ≈ JJ[F41,Z]
        @test saved_JJ[F41,Y]   ≈ JJ[F41,Y]
        @test saved_JJ[F41,G]   ≈ JJ[F41,G]
        @test saved_JJ[F41,TG]  ≈ JJ[F41,TG]
    end

    @testset "update lagged variables" begin
        @test saved_JJ[F19,LRRP] ≈ JJ[F19,LRRP]
        @test saved_JJ[F19,RR]   ≈ JJ[F19,RR]

        @test saved_JJ[F20,LIIP] ≈ JJ[F20,LIIP]
        @test saved_JJ[F20,II]   ≈ JJ[F20,II]

        @test saved_JJ[F24,LYP]  ≈ JJ[F24,LYP]
        @test saved_JJ[F24,Y]    ≈ JJ[F24,Y]

        @test saved_JJ[F31,LWP] ≈ JJ[F31,LWP]
        @test saved_JJ[F31,W]   ≈ JJ[F31,W]

        @test saved_JJ[F32,LXP] ≈ JJ[F32,LXP]
        @test saved_JJ[F32,X]   ≈ JJ[F32,X]
    end

    @testset "discount factor shock" begin
        @test saved_JJ[F33,BP] ≈ JJ[F33,BP]
        @test saved_JJ[F33,B]  ≈ JJ[F33,B]
    end

    @testset "g/y shock" begin
        @test saved_JJ[F34,GP] ≈ JJ[F34,GP]
        @test saved_JJ[F34,G]  ≈ JJ[F34,G]
    end

    @testset "tfp growth shock" begin
        @test saved_JJ[F35,ZP] ≈ JJ[F35,ZP]
        @test saved_JJ[F35,Z]  ≈ JJ[F35,Z]
    end

    @testset "investment shock" begin
        @test saved_JJ[F36,MUP] ≈ JJ[F36,MUP]
        @test saved_JJ[F36,MU]  ≈ JJ[F36,MU]
    end

    @testset "wage mkup shock" begin
        @test saved_JJ[F37,LAMWP] ≈ JJ[F37,LAMWP]
        @test saved_JJ[F37,LAMW]  ≈ JJ[F37,LAMW]
    end

    @testset "price mkup shock" begin
        @test saved_JJ[F38,LAMFP] ≈ JJ[F38,LAMFP]
        @test saved_JJ[F38,LAMF]  ≈ JJ[F38,LAMF]
    end

    @testset "monetary policy shock" begin
        @test saved_JJ[F39,MONP] ≈ JJ[F39,MONP]
        @test saved_JJ[F39,MON]  ≈ JJ[F39,MON]
    end
end

# Solution
if check_solution
    steadystate!(m, het)
    gx, hx, _ = klein(m)

    file = jldopen("$path/reference/solve.jld2", "r")
    saved_gx   = read(file, "gx")
    saved_hx   = read(file, "hx")
    close(file)

    @testset "Check solve outputs" begin
        @test saved_gx  ≈ gx
        @test saved_hx  ≈ hx
    end
end

# IRFs
if check_irfs
    endo = m.endogenous_states
    endo_aug = m.endogenous_states_augmented
    exo  = m.exogenous_shocks
    system = compute_system(m)
    states, obs, pseudo = impulse_responses(m, system, flip_shocks = true)

    file = jldopen("$path/reference/irfs.jld2", "r")
    IRFxp = read(file, "IRFxp")
    IRFC  = read(file, "IRFC")
    IRFyp = read(file, "IRFyp")

    nkp  = read(file, "nkp")
    ZP   = read(file, "ZP")
    ELLP = read(file, "ELLP")
    RRP  = read(file, "RRP")
    IIP  = read(file, "IIP")
    WP   = read(file, "WP")
    PIP  = read(file, "PIP")
    TTP  = read(file, "TTP")
    YP   = read(file, "YP")
    PIWP = read(file, "PIWP")
    XP   = read(file, "XP")
    KKP  = read(file, "KKP")
    QP   = read(file, "QP")
    RKP  = read(file, "RKP")
    HHP  = read(file, "HHP")
    TGP  = read(file, "TGP")
    BGP  = read(file, "BGP")
    MCP  = read(file, "MCP")
    LAMP = read(file, "LAMP")
    close(file)

    @testset "Check IRFs" begin
        @test states[endo[:z′_t], :, exo[:z_sh]] ≈ IRFxp[ZP, :]
        @test states[endo[:l′_t], :, exo[:z_sh]] ≈ IRFyp[ELLP - nkp, :]
        @test states[endo[:R′_t], :, exo[:z_sh]] ≈ IRFyp[RRP  - nkp, :]
        @test states[endo[:i′_t], :, exo[:z_sh]] ≈ IRFyp[IIP  - nkp, :]
        @test states[endo[:w′_t], :, exo[:z_sh]] ≈ IRFyp[WP   - nkp, :]
        @test states[endo[:π′_t], :, exo[:z_sh]] ≈ IRFyp[PIP  - nkp, :]
        @test states[endo[:t′_t], :, exo[:z_sh]] ≈ IRFyp[TTP  - nkp, :]
        @test states[endo[:y′_t], :, exo[:z_sh]] ≈ IRFyp[YP   - nkp, :]
        @test states[endo[:π_w′_t], :, exo[:z_sh]] ≈ IRFyp[PIWP - nkp, :]
        @test states[endo[:I′_t], :, exo[:z_sh]] ≈ IRFyp[XP   - nkp, :]
        @test states[endo[:k′_t], :, exo[:z_sh]] ≈ IRFxp[KKP, :]
        @test states[endo[:Q′_t], :, exo[:z_sh]] ≈ IRFyp[QP   - nkp, :]
        @test states[endo[:capreturn′_t], :, exo[:z_sh]] ≈ IRFyp[RKP - nkp, :]
        @test states[endo[:L′_t], :, exo[:z_sh]] ≈ IRFyp[HHP  - nkp, :]
        @test states[endo[:tg′_t], :, exo[:z_sh]] ≈ IRFyp[TGP - nkp, :]
        @test isapprox(states[endo[:bg′_t], :, exo[:z_sh]], IRFxp[BGP, :], atol = 1e-6)

        # Consumption IRF tests
        # contemporaneous consumption is the same as -l′_t
        @test states[endo[:l′_t], :, exo[:z_sh]] ≈ -IRFC
        @test states[endo_aug[:c_t1], :, exo[:z_sh]][2:end] ≈ IRFC[1:end-1]

        @test states[endo[:mc′_t], :, exo[:z_sh]] ≈ IRFyp[MCP - nkp, :]
        @test states[endo[:margutil′_t], :, exo[:z_sh]] ≈ IRFyp[LAMP - nkp, :]
    end
end
