using DSGE
using Test, BenchmarkTools
using JLD2

import DSGE: klein_transition_matrices, n_model_states, n_backward_looking_states

# What do you want to do?
check_steady_state = true
check_jacobian = true
check_solution = true

path = dirname(@__FILE__)

m = HetDSGE()

# Steady-state computation
if check_steady_state
    steadystate!(m)
    # @btime steadystate!(m)

    file = jldopen("$path/reference/steady_state.jld2", "r")
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

        @test saved_ell ≈ m[:lstar].value
        @test saved_c   ≈ m[:cstar].value
        @test saved_μ   ≈ m[:μstar].value
        # Tolerance of convergence of β is 1e-5
        @test abs(saved_β - m[:βstar].value) <= 1e-5
    end
end

if check_jacobian
    m.testing = true        # So that it will test against the unnormalized Jacobian
    JJ = DSGE.jacobian(m)
    # @btime JJ = DSGE.jacobian(m)

    file = jldopen("$path/reference/jacobian.jld2", "r")
    saved_JJ  = read(file, "JJ")
    LELLP = read(file, "LELLP")
    LMP = read(file, "LMP")
    KKP = read(file, "KKP")
    LRRP = read(file, "LRRP")
    LIIP = read(file, "LIIP")
    LPIP = read(file, "LPIP")
    L2PIP = read(file, "L2PIP")
    L3PIP = read(file, "L3PIP")
    LYP = read(file, "LYP")
    L2YP = read(file, "L2YP")
    L3YP = read(file, "L3YP")
    L4YP = read(file, "L4YP")
    LZP = read(file, "LZP")
    L2ZP = read(file, "L2ZP")
    L3ZP = read(file, "L3ZP")
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
    MP = read(file, "MP")
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
    LELL = read(file, "LELL")
    LM = read(file, "LM")
    KK = read(file, "KK")
    LRR = read(file, "LRR")
    LII = read(file, "LII")
    LPI = read(file, "LPI")
    L2PI = read(file, "L2PI")
    L3PI = read(file, "L3PI")
    LY = read(file, "LY")
    L2Y = read(file, "L2Y")
    L3Y = read(file, "L3Y")
    L4Y = read(file, "L4Y")
    LZ = read(file, "LZ")
    L2Z = read(file, "L2Z")
    L3Z = read(file, "L3Z")
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
    M = read(file, "M")
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
    F1 = read(file, "F1")
    F2 = read(file, "F2")
    F3 = read(file, "F3")
    F4 = read(file, "F4")
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
    F19 = read(file, "F19")
    F20 = read(file, "F20")
    F21 = read(file, "F21")
    F22 = read(file, "F22")
    F23 = read(file, "F23")
    F24 = read(file, "F24")
    F25 = read(file, "F25")
    F26 = read(file, "F26")
    F27 = read(file, "F27")
    F28 = read(file, "F28")
    F29 = read(file, "F29")
    F30 = read(file, "F30")
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

    endo = m.endogenous_states
    eq = m.equilibrium_conditions
    @testset "Check indices" begin
        @testset "function blocks which output a function" begin
            @test LELLP == endo[:l′_t1]  ## lagged ell function
            @test LMP   == endo[:μ′_t1]       ## lagged wealth distribution
        end
        @testset "endogenous scalar-valued states:" begin
            @test KKP   == first(endo[:k′_t])
            @test LRRP  == first(endo[:R′_t1])
            @test LIIP  == first(endo[:i′_t1])
            @test LPIP  == first(endo[:π′_t1])
            @test L2PIP == first(endo[:π′_t2])
            @test L3PIP == first(endo[:π′_t3])
            @test LYP   == first(endo[:y′_t1])
            @test L2YP  == first(endo[:y′_t2])
            @test L3YP  == first(endo[:y′_t3])
            @test L4YP  == first(endo[:y′_t4])
            @test LZP   == first(endo[:z′_t1])
            @test L2ZP  == first(endo[:z′_t2])
            @test L3ZP  == first(endo[:z′_t3])
            @test LWP   == first(endo[:w′_t1])
            @test LXP   == first(endo[:I′_t1])
        end
        @testset "# exogenous scalar-valued states" begin
            @test BP    == first(endo[:BP])
            @test GP    == first(endo[:GP])
            @test ZP    == first(endo[:ZP])
            @test MUP   == first(endo[:MUP])
            @test LAMWP == first(endo[:LAMWP])
            @test LAMFP == first(endo[:LAMFP])
            @test MONP  == first(endo[:MONP])
        end
        @testset "function-valued jumps" begin
            @test ELLP  == endo[:l′_t]
            @test MP    == endo[:μ′_t]
        end
        @testset "scalar-valued jumps" begin
            @test RRP   == first(endo[:R′_t])
            @test IIP   == first(endo[:i′_t])
            @test TTP   == first(endo[:t′_t])
            @test WP    == first(endo[:w′_t])
            @test HHP   == first(endo[:L′_t])
            @test PIP   == first(endo[:π′_t])
            @test PIWP  == first(endo[:wageinflation′_t])
            @test LAMP  == first(endo[:mu′_t])
            @test YP    == first(endo[:y′_t])
            @test XP    == first(endo[:I′_t])
            @test MCP   == first(endo[:mc′_t])
            @test QP    == first(endo[:Q′_t])
            @test RKP   == first(endo[:capreturn′_t])
            @test LELL  == endo[:l_t1]
            @test LM    == endo[:μ_t1]
            @test KK    == first(endo[:k_t])
            @test LRR   == first(endo[:R_t1])
            @test LII   == first(endo[:i_t1])
            @test LPI   == first(endo[:π_t1])
            @test L2PI  == first(endo[:π_t2])
            @test L3PI  == first(endo[:π_t3])
            @test LY    == first(endo[:y_t1])
            @test L2Y   == first(endo[:y_t2])
            @test L3Y   == first(endo[:y_t3])
            @test L4Y   == first(endo[:y_t4])
            @test LZ    == first(endo[:z_t1])
            @test L2Z   == first(endo[:z_t2])
            @test L3Z   == first(endo[:z_t3])
            @test LW    == first(endo[:w_t1])
            @test LX    == first(endo[:I_t1])
            @test B     == first(endo[:B])
            @test G     == first(endo[:G])
            @test Z     == first(endo[:Z])
            @test MU    == first(endo[:MU])
            @test LAMW  == first(endo[:LAMW])
            @test LAMF  == first(endo[:LAMF])
            @test MON   == first(endo[:MON])
            @test ELL   == endo[:ELL]
            @test M     == endo[:μ_t]
            @test RR    == first(endo[:R_t])
            @test II    == first(endo[:i_t])
            @test TT == first(endo[:t_t])
            @test W == first(endo[:w_t])
            @test HH == first(endo[:L_t])
            @test PI == first(endo[:π_t])
            @test PIW == first(endo[:wageinflation_t])
            @test LAM == first(endo[:mu_t])
            @test Y == first(endo[:y_t])
            @test X == first(endo[:I_t])
            @test MC == first(endo[:mc_t])
            @test Q == first(endo[:Q_t])
            @test RK == first(endo[:capreturn_t])
        end
        @testset "function blocks" begin
            @test F1  == eq[:eq_euler]
            @test F2  == eq[:eq_kolmogorov_fwd]
            @test F3  == eq[:eq_lag_ell]
            @test F4  == eq[:eq_lag_wealth]
        end
        @testset "function blocks which map functions to scalars" begin
            @test F5  == eq[:eq_market_clearing]
            @test F6  == eq[:eq_lambda]
        end
        @testset "# scalar blocks involving endogenous variables" begin
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
        @testset "lagged variables" begin
            @test F19 == eq[:LR]
            @test F20 == eq[:LI]
            @test F21 == eq[:LPI]
            @test F22 == eq[:L2PI]
            @test F23 == eq[:L3PI]
            @test F24 == eq[:LY]
            @test F25 == eq[:L2Y]
            @test F26 == eq[:L3Y]
            @test F27 == eq[:L4Y]
            @test F28 == eq[:LZ]
            @test F29 == eq[:L2Z]
            @test F30 == eq[:L3Z]
            @test F31 == eq[:LW]
            @test F32 == eq[:LX]
        end
        @testset "shocks" begin
            @test F33 == eq[:F33]
            @test F34 == eq[:F34]
            @test F35 == eq[:F35]
            @test F36 == eq[:F36]
            @test F37 == eq[:F37]
            @test F38 == eq[:F38]
            @test F39 == eq[:F39]
        end
    end

    @testset "Check jacobian outputs" begin
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
            @test saved_JJ[F2,LM] ≈ JJ[F2,LM]
            @test saved_JJ[F2,LRR] ≈ JJ[F2,LRR]
            @test saved_JJ[F2,LELL] ≈ JJ[F2,LELL]
            @test saved_JJ[F2,Z] ≈ JJ[F2,Z]
            @test saved_JJ[F2,M] ≈ JJ[F2,M]
            @test saved_JJ[F2,W] ≈ JJ[F2,W]
            @test saved_JJ[F2,HH] ≈ JJ[F2,HH]
            @test saved_JJ[F2,TT] ≈ JJ[F2,TT]
        end
        @testset "Check update lagged ELL" begin
            @test saved_JJ[F3,LELLP] ≈ JJ[F3,LELLP]
            @test saved_JJ[F3,ELL] ≈ JJ[F3,ELL]
        end
        @testset "Check update lagged M" begin
            @test saved_JJ[F4,LMP] ≈ JJ[F4,LMP]
            @test saved_JJ[F4,M] ≈ JJ[F4,M]
        end
        @testset "Check mkt clearing equation" begin
            @test saved_JJ[F5,Y] ≈ JJ[F5,Y]
            @test saved_JJ[F5,G] ≈ JJ[F5,G]
            @test saved_JJ[F5,X] ≈ JJ[F5,X]
            @test saved_JJ[F5,ELL] ≈ JJ[F5,ELL]
            @test saved_JJ[F5,M] ≈ JJ[F5,M]
        end
        @testset "Check lambda = average marginal utility" begin
            @test saved_JJ[F6,LAM] ≈ JJ[F6,LAM]
            @test saved_JJ[F6,M] ≈ JJ[F6,M]
            @test saved_JJ[F6,ELL] ≈ JJ[F6,ELL]
        end
        @testset "transfer equation" begin
            @test saved_JJ[F7,TT] ≈ JJ[F7,TT]
            @test saved_JJ[F7,RK] ≈ JJ[F7,RK]
            @test saved_JJ[F7,KK] ≈ JJ[F7,KK]
            @test saved_JJ[F7,Z] ≈ JJ[F7,Z]
            @test saved_JJ[F7,X] ≈ JJ[F7,X]
            @test saved_JJ[F7,MC] ≈ JJ[F7,MC]
            @test saved_JJ[F7,Y] ≈ JJ[F7,Y]
            @test saved_JJ[F7,G] ≈ JJ[F7,G]
        end
        @testset "investment" begin
            @test saved_JJ[F8,Q] ≈ JJ[F8,Q]
            @test saved_JJ[F8,MU] ≈ JJ[F8,MU]
            @test saved_JJ[F8,XP] ≈ JJ[F8,XP]
            @test saved_JJ[F8,ZP] ≈ JJ[F8,ZP]
            @test saved_JJ[F8,X] ≈ JJ[F8,X]
            @test saved_JJ[F8,LX] ≈ JJ[F8,LX]
            @test saved_JJ[F8,Z] ≈ JJ[F8,Z]
        end
        @testset "tobins q" begin
            @test saved_JJ[F9,RR] ≈ JJ[F9,RR]
            @test saved_JJ[F9,Q] ≈ JJ[F9,Q]
            @test saved_JJ[F9,RKP] ≈ JJ[F9,RKP]
            @test saved_JJ[F9,QP] ≈ JJ[F9,QP]
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
            @test saved_JJ[F16,LPI] ≈ JJ[F16,LPI]
            @test saved_JJ[F16,L2PI] ≈ JJ[F16,L2PI]
            @test saved_JJ[F16,L3PI] ≈ JJ[F16,L3PI]
            @test saved_JJ[F16,Y] ≈ JJ[F16,Y]
            @test saved_JJ[F16,L4Y] ≈ JJ[F16,L4Y]
            @test saved_JJ[F16,Z] ≈ JJ[F16,Z]
            @test saved_JJ[F16,LZ] ≈ JJ[F16,LZ]
            @test saved_JJ[F16,L2Z] ≈ JJ[F16,L2Z]
            @test saved_JJ[F16,L3Z] ≈ JJ[F16,L3Z]
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
        @testset "update lagged variables" begin
            @test saved_JJ[F19,LRRP] ≈ JJ[F19,LRRP]
            @test saved_JJ[F19,RR] ≈ JJ[F19,RR]
            @test saved_JJ[F20,LIIP] ≈ JJ[F20,LIIP]
            @test saved_JJ[F20,II] ≈ JJ[F20,II]
            @test saved_JJ[F21,LPIP] ≈ JJ[F21,LPIP]
            @test saved_JJ[F21,PI] ≈ JJ[F21,PI]
            @test saved_JJ[F22,L2PIP] ≈ JJ[F22,L2PIP]
            @test saved_JJ[F22,LPI] ≈ JJ[F22,LPI]
            @test saved_JJ[F23,L3PIP] ≈ JJ[F23,L3PIP]
            @test saved_JJ[F23,L2PI] ≈ JJ[F23,L2PI]
            @test saved_JJ[F24,LYP] ≈ JJ[F24,LYP]
            @test saved_JJ[F24,Y] ≈ JJ[F24,Y]
            @test saved_JJ[F25,L2YP] ≈ JJ[F25,L2YP]
            @test saved_JJ[F25,LY] ≈ JJ[F25,LY]
            @test saved_JJ[F26,L3YP] ≈ JJ[F26,L3YP]
            @test saved_JJ[F26,L2Y] ≈ JJ[F26,L2Y]
            @test saved_JJ[F27,L4YP] ≈ JJ[F27,L4YP]
            @test saved_JJ[F27,L3Y] ≈ JJ[F27,L3Y]
            @test saved_JJ[F28,LZP] ≈ JJ[F28,LZP]
            @test saved_JJ[F28,Z] ≈ JJ[F28,Z]
            @test saved_JJ[F29,L2ZP] ≈ JJ[F29,L2ZP]
            @test saved_JJ[F29,LZ] ≈ JJ[F29,LZ]
            @test saved_JJ[F30,L3ZP] ≈ JJ[F30,L3ZP]
            @test saved_JJ[F30,L2Z] ≈ JJ[F30,L2Z]
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
    m.testing = false
    gx, hx = klein(m)

    file = jldopen("$path/reference/solve.jld2", "r")
    saved_gx   = read(file, "gx")
    saved_hx   = read(file, "hx")
    close(file)

    @testset "Check solve outputs" begin
        @test saved_gx  ≈ gx
        @test saved_hx  ≈ hx
    end
end
