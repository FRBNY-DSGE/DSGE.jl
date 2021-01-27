using DSGE, Test, DelimitedFiles, JLD2, ModelConstructors, Random, FileIO
import DSGE: klein_transition_matrices, n_model_states, n_backward_looking_states

# What do you want to do?
check_steady_state = true
check_jacobian = true
check_solution = true
check_irfs = true
check_steady_state_calibrate = true
write_steady_state_calibrate = false
check_likelihood = true
write_likelihood = false

path = dirname(@__FILE__)

addtl_save_str = VERSION >= v"1.5" ? "_1.5" : ""

m = HetDSGEGovDebt(testing_gamma = true, ref_dir = HETDSGEGOVDEBT)
m <= Setting(:steady_state_only, true)
m <= Setting(:reduce_ell, true)
m <= Setting(:binsize, 10)

# Steady-state computation
if check_steady_state
    steadystate!(m)

    file = jldopen("$path/reference/steady_state_reduce_ell.jld2", "r")
    saved_Rk   = read(file, "Rk")
    saved_ω    = read(file, "omega")
    saved_kl   = read(file, "kl")
    saved_k    = read(file, "k")
    saved_y    = read(file, "y")
    saved_T    = read(file, "T")
    saved_ell  = read(file, "ell")
    saved_c    = read(file, "c")
    saved_μ    = read(file, "mu")
    saved_β    = read(file, "beta")
    close(file)

    @testset "Check steady state outputs" begin
        @test saved_Rk ≈ m[:Rkstar].value
        @test saved_ω  ≈ m[:ωstar].value
        @test saved_kl ≈ m[:klstar].value
        @test saved_k  ≈ m[:kstar].value
        @test saved_y  ≈ m[:ystar].value
        @test saved_T  ≈ m[:Tstar].value

        @test isapprox(saved_ell, m[:lstar].value, atol = .01)
        @test isapprox(saved_c, m[:cstar].value, atol = .01)
        @test isapprox(saved_μ, m[:μstar].value, atol = .01)
        # Tolerance of convergence of β is 1e-5
        @test isapprox(saved_β, m[:βstar].value, atol = .01)
    end
end

if check_jacobian
    m.testing = true # Ensures it tests against the unnormalized Jacobian
    m <= Setting(:steady_state_only, true)
    steadystate!(m)
    JJ = DSGE.jacobian(m)

    file = jldopen("$path/reference/jacobian_reduce_ell.jld2", "r")
    saved_JJ  = read(file, "JJ")
    KFP = read(file, "KFP")
    KKP = read(file, "KKP")
    LIIP = read(file, "LIIP")
    LYP = read(file, "LYP")
    LWP = read(file, "LWP")
    LXP = read(file, "LXP")
    BP = read(file, "BP")
    GP = read(file, "GP")
    ZP = read(file, "ZP")
    MUP = read(file, "MUP")
    LAMWP = read(file, "LAMWP")
    LAMFP = read(file, "LAMFP")
    MONP = read(file, "MONP")
    ELLP = read(file, "ELLP")
    CP = read(file, "CP")
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
    KF = read(file, "KF")
    KK = read(file, "KK")
    LII = read(file, "LII")
    LY = read(file, "LY")
    LW = read(file, "LW")
    LX = read(file, "LX")
    B = read(file, "B")
    G = read(file, "G")
    Z = read(file, "Z")
    MU = read(file, "MU")
    LAMW = read(file, "LAMW")
    LAMF = read(file, "LAMF")
    MON = read(file, "MON")
    ELL = read(file, "ELL")
    C = read(file, "C")
    RR = read(file, "RR")
    II = read(file, "II")
    TT = read(file, "TT")
    W = read(file, "W")
    HH = read(file, "HH")
    PI = read(file, "PI")
    PIW = read(file, "PIW")
    LAM = read(file, "LAM")
    Y = read(file, "Y")
    X = read(file, "X")
    MC = read(file, "MC")
    Q = read(file, "Q")
    RK = read(file, "RK")
    TG = read(file, "TG")
    TGP = read(file, "TGP")
    BG = read(file, "BG")
    BGP = read(file, "BGP")
    F1 = read(file, "F1")
    F2 = read(file, "F2")
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
    F40 = read(file, "F40")
    F41 = read(file, "F41")
    F42 = read(file, "F42")
    close(file)
    endo = DSGE.augment_model_states(m.endogenous_states_original,
                                     DSGE.n_model_states_original(m))

    eq = m.equilibrium_conditions
    @testset "Check indices" begin
        @testset "Function blocks which output a function" begin
            @test KFP == endo[:kf′_t]  ## lagged ell function
        end
        @testset "Endogenous scalar-valued states:" begin
            @test KKP   == first(endo[:k′_t])
            @test LIIP  == first(endo[:i′_t1])
            @test LYP   == first(endo[:y′_t1])
            @test LWP   == first(endo[:w′_t1])
            @test LXP   == first(endo[:I′_t1])
            @test BG == first(endo[:bg_t])
            @test BGP == first(endo[:bg′_t])
        end
        @testset "Exogenous scalar-valued states" begin
            @test BP    == first(endo[:b′_t])
            @test GP    == first(endo[:g′_t])
            @test ZP    == first(endo[:z′_t])
            @test MUP   == first(endo[:μ′_t])
            @test LAMWP == first(endo[:λ_w′_t])
            @test LAMFP == first(endo[:λ_f′_t])
            @test MONP  == first(endo[:rm′_t])
        end
        @testset "Function-valued jumps" begin
            @test ELLP  == endo[:l′_t]
        end
        @testset "Scalar-valued jumps" begin
            @test CP   == first(endo[:C′_t])
            @test RRP  == first(endo[:R′_t])
            @test IIP  == first(endo[:i′_t])
            @test TTP  == first(endo[:t′_t])
            @test WP   == first(endo[:w′_t])
            @test HHP  == first(endo[:L′_t])
            @test PIP  == first(endo[:π′_t])
            @test PIWP == first(endo[:π_w′_t])
            @test LAMP == first(endo[:margutil′_t])
            @test YP   == first(endo[:y′_t])
            @test XP   == first(endo[:I′_t])
            @test MCP  == first(endo[:mc′_t])
            @test QP   == first(endo[:Q′_t])
            @test RKP  == first(endo[:capreturn′_t])
            @test KF   == endo[:kf_t]
            @test KK   == first(endo[:k_t])
            @test LII  == first(endo[:i_t1])
            @test LY   == first(endo[:y_t1])
            @test LW   == first(endo[:w_t1])
            @test LX   == first(endo[:I_t1])
            @test B    == first(endo[:b_t])
            @test G    == first(endo[:g_t])
            @test Z    == first(endo[:z_t])
            @test MU   == first(endo[:μ_t])
            @test LAMW == first(endo[:λ_w_t])
            @test LAMF == first(endo[:λ_f_t])
            @test MON  == first(endo[:rm_t])
            @test ELL  == endo[:l_t]
            @test C    == first(endo[:C_t])
            @test RR   == first(endo[:R_t])
            @test II   == first(endo[:i_t])
            @test TT   == first(endo[:t_t])
            @test W    == first(endo[:w_t])
            @test HH   == first(endo[:L_t])
            @test PI   == first(endo[:π_t])
            @test PIW  == first(endo[:π_w_t])
            @test LAM  == first(endo[:margutil_t])
            @test Y    == first(endo[:y_t])
            @test X    == first(endo[:I_t])
            @test MC   == first(endo[:mc_t])
            @test Q    == first(endo[:Q_t])
            @test RK   == first(endo[:capreturn_t])
            @test TG   == first(endo[:tg_t])
            @test TGP  == first(endo[:tg′_t])
        end
        @testset "Function blocks" begin
            @test F1  == eq[:eq_euler]
            @test F2  == eq[:eq_kolmogorov_fwd]
        end
        @testset "Function blocks which map functions to scalars" begin
            @test F5  == eq[:eq_agg_consumption]
            @test F6  == eq[:eq_lambda]
        end
        @testset "Scalar blocks involving endogenous variables" begin
            @test F7  == eq[:eq_transfers]
            @test F8  == eq[:eq_investment]
            @test F9  == eq[:eq_tobin_q]
            @test F10 == eq[:eq_capital_accumulation]
            @test F11 == eq[:eq_wage_phillips]
            @test F12 == eq[:eq_price_phillips]
            @test F13 == eq[:eq_marginal_cost]
            @test F14 == eq[:eq_gdp]
            @test F15 == eq[:eq_optimal_kl]
            @test F16 == eq[:eq_taylor]
            @test F17 == eq[:eq_fisher]
            @test F18 == eq[:eq_nominal_wage_inflation]
        end
        @testset "Lagged variables" begin
            @test F20 == eq[:LI]
            @test F24 == eq[:LY]
            @test F31 == eq[:LW]
            @test F32 == eq[:LX]
        end
        @testset "Shocks" begin
            @test F33 == eq[:eq_b]
            @test F34 == eq[:eq_g]
            @test F35 == eq[:eq_z]
            @test F36 == eq[:eq_μ]
            @test F37 == eq[:eq_λ_w]
            @test F38 == eq[:eq_λ_f]
            @test F39 == eq[:eq_rm]
        end
        @testset "New equations" begin
            @test F40 == eq[:eq_fiscal_rule]
            @test F41 == eq[:eq_g_budget_constraint]
            @test F42 == eq[:eq_resource_constraint]
        end
    end

    @testset "Check Jacobian outputs" begin
        @testset "Euler Equation" begin
            @test saved_JJ[F1,ELLP] ≈ JJ[F1, ELLP]
            @test saved_JJ[F1,ZP] ≈ JJ[F1,ZP]
            @test saved_JJ[F1,WP] ≈ JJ[F1,WP]
            @test saved_JJ[F1,HHP] ≈ JJ[F1,HHP]
            @test saved_JJ[F1,TTP] ≈ JJ[F1,TTP]
            @test saved_JJ[F1,B] ≈ JJ[F1,B]
            @test saved_JJ[F1,ELL] ≈ JJ[F1,ELL]
            @test saved_JJ[F1,RR] ≈ JJ[F1,RR]
        end
        @testset "Check Kolmogorov Fwd Equation" begin
            @test saved_JJ[F2,KFP] ≈ JJ[F2,KFP]
            @test saved_JJ[F2,KF] ≈ JJ[F2,KF]
            @test saved_JJ[F2,ELL] ≈ JJ[F2,ELL]
            @test saved_JJ[F2,RR] ≈ JJ[F2,RR]
            @test saved_JJ[F2,Z] ≈ JJ[F2,Z]
            @test saved_JJ[F2,W] ≈ JJ[F2,W]
            @test saved_JJ[F2,HH] ≈ JJ[F2,HH]
            @test saved_JJ[F2,TT] ≈ JJ[F2,TT]
        end
        @testset "Check aggregate consumption" begin
            @test saved_JJ[F5,C]   ≈ JJ[F5,C]
            @test saved_JJ[F5,ELL] ≈ JJ[F5,ELL]
            @test saved_JJ[F5,KF]  ≈ JJ[F5,KF]
            @test saved_JJ[F5,Z]   ≈ JJ[F5,Z]
            @test saved_JJ[F5,W]   ≈ JJ[F5,W]
            @test saved_JJ[F5,HH]  ≈ JJ[F5,HH]
            @test saved_JJ[F5,TT]  ≈ JJ[F5,TT]
        end
        @testset "Check lambda = average marginal utility" begin
            @test saved_JJ[F6,LAM] ≈ JJ[F6,LAM]
            @test saved_JJ[F6,KF]  ≈ JJ[F6,KF]
            @test saved_JJ[F6,RR]  ≈ JJ[F6,RR]
            @test saved_JJ[F6,Z]   ≈ JJ[F6,Z]
            @test saved_JJ[F6,W]   ≈ JJ[F6,W]
            @test saved_JJ[F6,HH]  ≈ JJ[F6,HH]
            @test saved_JJ[F6,TT]  ≈ JJ[F6,TT]
            @test saved_JJ[F6,ELL] ≈ JJ[F6,ELL]
        end
        @testset "transfer equation" begin
            @test saved_JJ[F7,TT] ≈ JJ[F7,TT]
            @test saved_JJ[F7,RK] ≈ JJ[F7,RK]
            @test saved_JJ[F7,KK] ≈ JJ[F7,KK]
            @test saved_JJ[F7,Z]  ≈ JJ[F7,Z]
            @test saved_JJ[F7,X]  ≈ JJ[F7,X]
            @test saved_JJ[F7,MC] ≈ JJ[F7,MC]
            @test saved_JJ[F7,Y]  ≈ JJ[F7,Y]
            @test saved_JJ[F7,G]  ≈ JJ[F7,G]
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
        @testset "tobins q" begin
            @test saved_JJ[F9,RR]  ≈ JJ[F9,RR]
            @test saved_JJ[F9,Q]   ≈ JJ[F9,Q]
            @test saved_JJ[F9,RKP] ≈ JJ[F9,RKP]
            @test saved_JJ[F9,QP]  ≈ JJ[F9,QP]
        end
        @testset "capital accumulation" begin
            @test saved_JJ[F10,KKP] ≈ JJ[F10,KKP]
            @test saved_JJ[F10,KK] ≈ JJ[F10,KK]
            @test saved_JJ[F10,Z] ≈ JJ[F10,Z]
            @test saved_JJ[F10,MU] ≈ JJ[F10,MU]
            @test saved_JJ[F10,X] ≈ JJ[F10,X]
        end
        @testset "wage phillips curve" begin
            @test saved_JJ[F11,PIW] ≈ JJ[F11,PIW]
            @test saved_JJ[F11,LAMW] ≈ JJ[F11,LAMW]
            @test saved_JJ[F11,HH] ≈ JJ[F11,HH]
            @test saved_JJ[F11,LAM] ≈ JJ[F11,LAM]
            @test saved_JJ[F11,W] ≈ JJ[F11,W]
            @test saved_JJ[F11,PIP] ≈ JJ[F11,PIP]
        end
        @testset "price phillips curve" begin
            @test saved_JJ[F12,PI] ≈ JJ[F12,PI]
            @test saved_JJ[F12,MC] ≈ JJ[F12,MC]
            @test saved_JJ[F12,LAMF] ≈ JJ[F12,LAMF]
            @test saved_JJ[F12,PIP] ≈ JJ[F12,PIP]
        end
        @testset "marginal cost" begin
            @test saved_JJ[F13,MC] ≈ JJ[F13,MC]
            @test saved_JJ[F13,W] ≈ JJ[F13,W]
            @test saved_JJ[F13,RK] ≈ JJ[F13,RK]
        end
        @testset "gdp" begin
            @test saved_JJ[F14,Y] ≈ JJ[F14,Y]
            @test saved_JJ[F14,Z] ≈ JJ[F14,Z]
            @test saved_JJ[F14,KK] ≈ JJ[F14,KK]
            @test saved_JJ[F14,HH] ≈ JJ[F14,HH]
        end
        @testset "optimal k/l ratio" begin
            @test saved_JJ[F15,RK] ≈ JJ[F15,RK]
            @test saved_JJ[F15,W] ≈ JJ[F15,W]
            @test saved_JJ[F15,HH] ≈ JJ[F15,HH]
            @test saved_JJ[F15,KK] ≈ JJ[F15,KK]
            @test saved_JJ[F15,Z] ≈ JJ[F15,Z]
        end
        @testset "taylor rule" begin
            @test saved_JJ[F16,II] ≈ JJ[F16,II]
            @test saved_JJ[F16,LII] ≈ JJ[F16,LII]
            @test saved_JJ[F16,PI] ≈ JJ[F16,PI]
            @test saved_JJ[F16,Y] ≈ JJ[F16,Y]
            @test saved_JJ[F16,LY] ≈ JJ[F16,LY]
            @test saved_JJ[F16,Z] ≈ JJ[F16,Z]
            @test saved_JJ[F16,MON] ≈ JJ[F16,MON]
        end
        @testset "fisher eqn" begin
            @test saved_JJ[F17,RR] ≈ JJ[F17,RR]
            @test saved_JJ[F17,PIP] ≈ JJ[F17,PIP]
            @test saved_JJ[F17,II] ≈ JJ[F17,II]
        end
       @testset "wage inflation" begin
           @test saved_JJ[F18,PIW] ≈ JJ[F18,PIW]
           @test saved_JJ[F18,PI] ≈ JJ[F18,PI]
           @test saved_JJ[F18,Z] ≈ JJ[F18,Z]
           @test saved_JJ[F18,W] ≈ JJ[F18,W]
           @test saved_JJ[F18,LW] ≈ JJ[F18,LW]
       end
       @testset "fiscal rule" begin
           @test saved_JJ[F40,TG] ≈ JJ[F40,TG]
           @test saved_JJ[F40,RR] ≈ JJ[F40,RR]
           @test saved_JJ[F40,BG] ≈ JJ[F40,BG]
           @test saved_JJ[F40,Z] ≈ JJ[F40,Z]
           @test saved_JJ[F40,Y] ≈ JJ[F40,Y]
           @test saved_JJ[F40,G] ≈ JJ[F40,G]
       end
       @testset "govt budget constraint" begin
           @test saved_JJ[F41,BGP] ≈ JJ[F41,BGP]
           @test saved_JJ[F41,RR] ≈ JJ[F41,RR]
           @test saved_JJ[F41,BG] ≈ JJ[F41,BG]
           @test saved_JJ[F41,Z] ≈ JJ[F41,Z]
           @test saved_JJ[F41,Y] ≈ JJ[F41,Y]
           @test saved_JJ[F41,G] ≈ JJ[F41,G]
           @test saved_JJ[F41,TG] ≈ JJ[F41,TG]
       end
       @testset "update lagged variables" begin
            @test saved_JJ[F20,LIIP] ≈ JJ[F20,LIIP]
            @test saved_JJ[F20,II] ≈ JJ[F20,II]
            @test saved_JJ[F24,LYP] ≈ JJ[F24,LYP]
            @test saved_JJ[F24,Y] ≈ JJ[F24,Y]
            @test saved_JJ[F31,LWP] ≈ JJ[F31,LWP]
            @test saved_JJ[F31,W] ≈ JJ[F31,W]
            @test saved_JJ[F32,LXP] ≈ JJ[F32,LXP]
            @test saved_JJ[F32,X] ≈ JJ[F32,X]
        end
        @testset "discount factor shock" begin
            @test saved_JJ[F33,BP] ≈ JJ[F33,BP]
            @test saved_JJ[F33,B] ≈ JJ[F33,B]
        end
        @testset "g/y shock" begin
            @test saved_JJ[F34,GP] ≈ JJ[F34,GP]
            @test saved_JJ[F34,G] ≈ JJ[F34,G]
        end
        @testset " tfp growth shock" begin
            @test saved_JJ[F35,ZP] ≈ JJ[F35,ZP]
            @test saved_JJ[F35,Z] ≈ JJ[F35,Z]
        end
        @testset "investment shock" begin
            @test saved_JJ[F36,MUP] ≈ JJ[F36,MUP]
            @test saved_JJ[F36,MU] ≈ JJ[F36,MU]
        end
        @testset "wage mkup shock" begin
            @test saved_JJ[F37,LAMWP] ≈ JJ[F37,LAMWP]
            @test saved_JJ[F37,LAMW] ≈ JJ[F37,LAMW]
        end
        @testset "price mkup shock" begin
            @test saved_JJ[F38,LAMFP] ≈ JJ[F38,LAMFP]
            @test saved_JJ[F38,LAMF] ≈ JJ[F38,LAMF]
        end
        @testset "monetary policy shock" begin
            @test saved_JJ[F39,MONP] ≈ JJ[F39,MONP]
            @test saved_JJ[F39,MON] ≈ JJ[F39,MON]
        end
    end
end

if check_solution
    if !check_steady_state
        steadystate!(m)
    end
    nx_save = get_setting(m, :nx)
    m.testing = false
    m <= Setting(:nx, nx_save)
    gx, hx, _ = klein(m)

    file = JLD2.jldopen("$path/reference/solve_reduce_ell.jld2", "r")
    saved_gx   = read(file, "gx")
    saved_hx   = read(file, "hx")
    close(file)

    @testset "Check solve outputs" begin
        @test saved_gx  ≈ gx
        @test saved_hx  ≈ hx
    end
end

if check_irfs
    steadystate!(m)
    file = jldopen("$path/reference/irfs_reduce_ell.jld2", "r")

    IRFkf = read(file, "IRFkf")
    IRFZ = read(file, "IRFZ")
    IRFsh = read(file, "IRFsh")
    IRFell = read(file, "IRFell")
    IRFR = read(file, "IRFR")
    IRFi = read(file, "IRFi")
    IRFW = read(file, "IRFW")
    IRFPI = read(file, "IRFPI")
    IRFT = read(file, "IRFT")
    IRFY = read(file, "IRFY")
    IRFPIW = read(file, "IRFPIW")
    IRFX = read(file, "IRFX")
    IRFK = read(file, "IRFK")
    IRFQ = read(file, "IRFQ")
    IRFRK = read(file, "IRFRK")
    IRFH = read(file, "IRFH")
    IRFc = read(file, "IRFc")
    IRFC = read(file, "IRFC")
    IRFMC = read(file, "IRFMC")
    IRFLAM = read(file, "IRFLAM")
    IRFm = read(file, "IRFm")
    close(file)

    sys = compute_system(m)
    states, obs, pseudo = impulse_responses(m, sys, use_augmented_states = true)
    endo = m.endogenous_states_original

    @testset "Check IRFs" begin
        #Last entry is 3 because it's the third shock, the Z shock
        @test IRFkf ≈ states[endo[:kf′_t], 1:20, 3]
        @test IRFsh ≈ vec(states[endo[:z′_t], 1:20, 3])
        @test IRFell ≈ states[endo[:l′_t], 1:20, 3]
        @test IRFR ≈ vec(states[endo[:R′_t], 1:20, 3])
        @test IRFi ≈ vec(states[endo[:i′_t], 1:20, 3])
        @test IRFW ≈ vec(states[endo[:w′_t], 1:20, 3])
        @test IRFPI ≈ vec(states[endo[:π′_t], 1:20, 3])
        @test IRFT ≈ vec(states[endo[:t′_t], 1:20, 3])
        @test IRFY ≈ vec(states[endo[:y′_t], 1:20, 3])
        @test IRFPIW ≈ vec(states[endo[:π_w′_t], 1:20, 3])
        @test IRFX ≈ vec(states[endo[:I′_t], 1:20, 3])
        @test IRFK ≈ vec(states[endo[:k′_t], 1:20, 3])
        @test IRFQ ≈ vec(states[endo[:Q′_t], 1:20, 3])
        @test IRFRK ≈ vec(states[endo[:capreturn′_t], 1:20, 3])
        @test IRFH ≈ vec(states[endo[:L′_t], 1:20, 3])
        #@test IRFC ≈ vec(obs[m.observables[:obs_consumption], 1:20, 3])
        @test IRFMC ≈ vec(states[endo[:mc′_t], 1:20, 3])
        @test IRFLAM ≈ vec(states[endo[:margutil′_t], 1:20, 3])
        #@test IRFm ≈ vec(states[endo[:margutil′_t], 1:20, 3])
    end
end

m <= Setting(:steady_state_only, false)

# Steady-state computation
if check_steady_state_calibrate || write_steady_state_calibrate
    Random.seed!(0)
    steadystate!(m)

    if write_steady_state_calibrate
        JLD2.jldopen("$path/reference/steady_state_calibration$(addtl_save_str).jld2", "w") do file
            file["c"]    = m[:cstar].value
            file["m"]    = m[:μstar].value
            file["mpc"]  = m[:mpc].value
            file["pc0"]  = m[:pc0].value
            file["beta"] = m[:βstar].value
        end
    end

    if check_steady_state_calibrate
        file = jldopen("$path/reference/steady_state_calibration$(addtl_save_str).jld2", "r")
        saved_c    = read(file, "c")
        saved_μ    = read(file, "m")
        saved_mpc  = read(file, "mpc")
        saved_β    = read(file, "beta")
        saved_pc0  = read(file, "pc0")
        close(file)

        @testset "Check steady state outputs" begin
            @test saved_c   ≈ m[:cstar].value
            @test saved_μ   ≈ m[:μstar].value
            @test saved_pc0 ≈ m[:pc0].value
            @test saved_mpc ≈ m[:mpc].value
            # Tolerance of convergence of β is 1e-5
            @test saved_β ≈ m[:βstar].value
        end
    end
end

if check_likelihood || write_likelihood
    data = Matrix{Float64}(readdlm("$path/../../../reference/YY.txt")')

    if write_likelihood
        println("Overwriting likelihood...")
        JLD2.jldopen("$path/reference/likelihood_reduce_ell$(addtl_save_str).jld2", true, true,
                     true, IOStream) do file
            file["likelihood"] = likelihood(m, data)
        end
    end
    if check_likelihood
        lik_save = load("$path/reference/likelihood_reduce_ell$(addtl_save_str).jld2", "likelihood")

        @testset "Checking Likelihood output" begin
            @test likelihood(m, data) ≈ lik_save
        end
    end
end
