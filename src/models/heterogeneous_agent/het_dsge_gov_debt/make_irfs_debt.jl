#close("all")
dur =20 #Periods to compute IRFs

(nkc, nkp) = size(Qx) # nkc = number of predetermined variables in  `coefficient values'
                      # representation (i.e. after removing degrees of freedom to normalize
                      # the distribution - this is what we use to solve the model etc.)
                      # nkp = number of predetermined variables in `grid values' representation
                      # (i.e. values on a grid, which is what we use to make pictures)
df = nkp-nkc # missing `degrees of freedom'

whichshock = ZP

if whichshock==BP
	shname = "b shock"
elseif whichshock==GP
	shname = "govt spending"
elseif whichshock==ZP
	shname = "TFP growth"
elseif whichshock==MUP
	shname = "MEI shock"
elseif whichshock==LAMFP
	shname = "markup shock"
elseif whichshock==LAMWP
	shname = "wage markup"
elseif whichshock==MONP
	shname = "monetary pol"
end

shock        = -0.15
IRFxc        = zeros(nkc,dur) #IRF to income shock, coefficient values
IRFxc[whichshock - df,1] = shock # ZP - df is the location of TFP shock in coefficient values rep
for t=1:dur-1
    IRFxc[:,t+1] = hx*IRFxc[:,t]
end
IRFyc = gx*IRFxc #y response

#Shocks in terms of grid values
IRFxp  = Qx'IRFxc
IRFyp = Qy'IRFyc

IRFkf   = IRFxp[KFP,:]
IRFZ   = IRFxp[ZP,:]
IRFsh  = IRFxp[whichshock,:]
IRFell = IRFyp[ELLP - nkp,:] # ELLP - nkp is the location of ELLP in zIRFyp
IRFR   = IRFyp[RRP - nkp,:]
IRFi   = IRFyp[IIP - nkp,:]
IRFW   = IRFyp[WP - nkp,:]
IRFPI  = IRFyp[PIP - nkp,:]
IRFT   = IRFyp[TTP - nkp,:]
IRFY   = IRFyp[YP - nkp,:]
IRFPIW   = IRFyp[PIWP - nkp,:]
IRFX   = IRFyp[XP - nkp,:]
IRFK   = IRFxp[KKP,:]
IRFQ   = IRFyp[QP - nkp,:]
IRFRK   = IRFyp[RKP - nkp,:]
IRFH   = IRFyp[HHP - nkp,:]
IRFTG   = IRFyp[TGP - nkp,:]
IRFBG = IRFxp[BGP,:]
IRFm = IRFkf + dF2_dRZ*(IRFR' - IRFZ') + dF2_dWH*(IRFW'+IRFH') + dF2_dTT*IRFT'


IRFc   = -repmat(unc.*c,1,dur).*IRFell
IRFC   = ((m.*aswts)'*IRFc)' + ((aswts.*c)'*IRFm)' # in levels
IRFMC  = IRFyp[MCP - nkp,:]
IRFLAM  = IRFyp[LAMP - nkp,:]

if save_data
    file = JLD2.jldopen("irfs.jld2", "w")
    write(file, "IRFkf", IRFkf)
    write(file, "IRFZ", IRFZ)
    write(file, "IRFsh", IRFsh)
    write(file, "IRFell", IRFell)
    write(file, "IRFR", IRFR)
    write(file, "IRFi", IRFi)
    write(file, "IRFW", IRFW)
    write(file, "IRFPI", IRFPI)
    write(file, "IRFT", IRFT)
    write(file, "IRFY", IRFY)
    write(file, "IRFPIW", IRFPIW)
    write(file, "IRFX", IRFX)
    write(file, "IRFK", IRFK)
    write(file, "IRFQ", IRFQ)
    write(file, "IRFRK", IRFRK)
    write(file, "IRFH", IRFH)
    write(file, "IRFc", IRFc)
    write(file, "IRFC", IRFC)
    write(file, "IRFMC", IRFMC)
    write(file, "IRFLAM", IRFLAM)
    write(file, "IRFm", IRFm)
    write(file, "IRFTG", IRFTG)
    write(file, "IRFBG", IRFBG)
    close(file)
end

#=
figure()
plot(agrid,c[1:na],agrid,c[na+1:2*na])

figure()
plot(agrid,m[1:na],agrid,m[na+1:2*na])

fig = figure(figsize=(12,6))

subplot(3,6,1)
grid("on")
plot(1:dur,IRFsh)
title(shname)

subplot(3,6,2)
grid("on")
plot(1:dur,IRFR)
title("Real interest rate")

subplot(3,6,3)
grid("on")
plot(1:dur,IRFi)
title("Nominal interest rate")

subplot(3,6,4)
grid("on")
plot(1:dur,IRFPI)
title("Inflation")

subplot(3,6,5)
grid("on")
plot(1:dur,IRFC)
title("Consumption")

subplot(3,6,6)
grid("on")
plot(1:dur,IRFLAM)
title("Lambda")

subplot(3,6,7)
grid("on")
plot(1:dur,IRFW)
title("Real wage")

subplot(3,6,8)
grid("on")
plot(1:dur,IRFPIW)
title("Nominal Wage Inflation")



subplot(3,6,9)
grid("on")
plot(1:dur,IRFT)
title("Transfers")

subplot(3,6,10)
grid("on")
plot(1:dur,IRFQ)
title("Tobin's Q")

subplot(3,6,11)
grid("on")
plot(1:dur,IRFMC)
title("Marginal cost")

subplot(3,6,12)
grid("on")
plot(1:dur,IRFRK)
title("Return on capital")





subplot(3,6,13)
grid("on")
plot(1:dur,IRFX)
title("Investment")

subplot(3,6,14)
grid("on")
plot(1:dur,IRFY)
title("Output")

subplot(3,6,15)
grid("on")
plot(1:dur,IRFK)
title("Capital")

subplot(3,6,16)
grid("on")
plot(1:dur,IRFH)
title("Hours")

subplot(3,6,17)
grid("on")
plot(1:dur,IRFBG)
title("Govt Debt")

subplot(3,6,18)
grid("on")
plot(1:dur,IRFTG)
title("Taxes")



# distributional irfs

ag = agrid # might truncate; for now we don't

thingtoplot = IRFm[1:na,:];
fig = figure(figsize=(8,6))
#ax = fig[:gca](projection="3d")
ax = Axes3D(fig)
ax[:set_zlim](minimum(thingtoplot), maximum(thingtoplot));
xgridgph = repmat(ag,1,dur);
ygridgph = repmat((1:dur)',na,1);
ax[:plot_surface](xgridgph, ygridgph, thingtoplot, cmap=ColorMap("jet"),
                  edgecolors="k",linewidth=0.25)
xlabel("a")
ylabel("t")
title("IRF of m, low skill")

thingtoplot = IRFm[na+1:2*na,:];
fig = figure(figsize=(8,6))
#ax = fig[:gca](projection="3d")
ax = Axes3D(fig)
ax[:set_zlim](minimum(thingtoplot), maximum(thingtoplot));
xgridgph = repmat(ag,1,dur);
ygridgph = repmat((1:dur)',na,1);
ax[:plot_surface](xgridgph, ygridgph, thingtoplot, cmap=ColorMap("jet"),
                  edgecolors="k",linewidth=0.25)
xlabel("a")
ylabel("t")
title("IRF of m, high skill")

thingtoplot = IRFc[1:na,:];
fig = figure(figsize=(8,6))
#ax = fig[:gca](projection="3d")
ax = Axes3D(fig)
ax[:set_zlim](minimum(thingtoplot), maximum(thingtoplot));
xgridgph = repmat(ag,1,dur);
ygridgph = repmat((1:dur)',na,1);
ax[:plot_surface](xgridgph, ygridgph, thingtoplot, cmap=ColorMap("jet"),
                  edgecolors="k",linewidth=0.25)
xlabel("a")
ylabel("t")
title("IRF of c, low skill")

thingtoplot = IRFc[na+1:2*na,:];
fig = figure(figsize=(8,6))
#ax = fig[:gca](projection="3d")
ax = Axes3D(fig)
ax[:set_zlim](minimum(thingtoplot), maximum(thingtoplot));
xgridgph = repmat(ag,1,dur);
ygridgph = repmat((1:dur)',na,1);
ax[:plot_surface](xgridgph, ygridgph, thingtoplot, cmap=ColorMap("jet"),
                  edgecolors="k",linewidth=0.25)
xlabel("a")
ylabel("t")
title("IRF of c, high skill")

=#
