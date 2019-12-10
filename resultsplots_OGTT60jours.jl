# NEEDED MODULES  -- This can take time at the first launch
# Code main be needed to run twice to plot correctly (shift-enter to relaunch the cell)
# The second time is much faster
using DifferentialEquations
using Plots
using LaTeXStrings
# Load the backend (engine) for plotting based on javascript
plotlyjs() # FOR INTERACTIVE PLOT
# plotly()   

function weighted_leastsquare(O,R,E)
    denom  = sum(O.^2)
    weights = 1. / (E.^2)
    return sqrt(sum(weights.*(O-R).^2)/denom)
end

function leastsquare(O,R)
    denom  = sum(O.^2)
    return sqrt(sum((O-R).^2)/denom)
end

# Saturation function to positive values 
# (concentrations are positive function)
# This is just used as a safety for the parameter search phase
@inline function sat(x)
    if x < 0.0
        return 0.0
    else
        x
    end
end

# Gaussian model of plasma Glucose appearance rate
function GlcGauss(t)
    σ = 35.0 # width
    μ = 30.0 # peak time
    dose = 3982 # total glucose dose (area under curve)
    return dose*1/(σ*sqrt(2*pi))*exp((-(t-μ)^2/(2*σ^2)))
end

# Parametric ODEs of te Glucose response model
# t : time
# u : variables
# p : parameters
# Here I call the gaussian model of glucose input GlcGauss(t)
# And I saturate the variables to the positive cone
function f_param!(du,u,p,t)
 du[1] = -p[1]*(sat(u[1])-p[2])- p[10]*sat(u[2])*sat(u[1])+GlcGauss(t)
 du[2] = -p[3]*sat(u[2])+p[4]*(sat(u[3])-p[5])
 du[3] = -p[6]*sat(u[3])+p[7]*(sat(u[1])-p[8])*t + p[9]*u[4]
 du[4] = -p[9]*u[4]           
end   
  

# Initial conditions at 60j
u0 = [95.0,0.0,34.0,5950.0]
tspan = (0.0,120.0)


# #  A Parameter set from the MINMOD model as a reference
# #       p1  Gb    p2    p3     Ib  n   γ    h     |added parameters => |     p4    r
# p0 = [0.01,100.0,0.01,0.00001,8.0,0.1,0.01,100.0,                           0.07,1.0]

    

# Fitted parameter set to the control group at 60 Days ===> err: 0.00251
#        p1    Gb     p2     p3    Ib     n       γ     h     p4      r
  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]    
###################################################
## Simulation of the Control experiment at 26Days      
prob_param = ODEProblem(f_param!,u0,tspan,p0)  
sctrl = solve(prob_param, alg_hints = [:stiff])
## Experimental results    
texp = [0,10,20,30,45,60,90,120]
glcexp = [91.9,155.5,186.0,204.9,197.2,168.7,130.5,108.2]
uppererror = [2.9,9.1,7.0,6.6,6.5,6.7,2.9,2.9]
lowererror = [2.9,9.1,7.0,6.6,6.5,6.7,2.9,2.9]
errs = [lowererror;uppererror]
res_tot = [sctrl(ti) for ti in texp]
res = [ri[1] for ri in res_tot]
## Computation of the weighted least square error and least square
wls_err_Ctrl = weighted_leastsquare(glcexp,res,uppererror)
print("CONTROL wls = ", wls_err_Ctrl, "\n")
ls_err_Ctrl = leastsquare(glcexp,res)
## Generation of the figure and first plots 
plot(sctrl,vars=[(0,1)],label="Glucose (control) simulations")
scatter!(texp,glcexp,yerr = errs,marker = stroke(2, :blue),
       xlabel="time (min)", ylabel= "Concentration: mg/dL",
    label="Control exp results", legend=:topright, markerstrokecolor = :auto)

    
    
# Parameters sets for the different hypothesis for the group CD1 at 60 Days
########################################################################    
# Control parameter set  ===>   err 0.00440
#        p1    Gb     p2     p3    Ib     n       γ      h     p4      r
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 1.1: sensibility to interstitial insuline (r) ===> err : 0.00395
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  0.97]
# hypothese 1.2: disponibility of insuline for X(t) (p3) ==> err: 0.00395
#  p0 = [0.01, 100.0, 0.79, 0.0325, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 2: degradation  insuline (n) ==> err : 0.00374
# p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.49, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 3: interstitial insuline degradation (p2) ==> err : 0.00396
#  p0 = [0.01, 100.0, 0.82, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothesis 4: γ (Insuline production acceleration) decreased  ==> err =  0.00267
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0072, 65.0, 0.0170,  1.0]
# hypothesis 5:  change p4 (stored insuline release rate)  ==> err = 0.00439       
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0168,  1.0]
# hypothese 6: reponse Glucose threshold (h)  ==> err:  0.00204
  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 73.0, 0.0170,  1.0]
########################################################################
# Simulation CD1 60 days    
prob_param = ODEProblem(f_param!,u0,tspan,p0)
scd1 = solve(prob_param, alg_hints = [:stiff])
# Experimental results    
glcexpCd1 = [96.4,154.9,190.9,207.8,199.9,173.2,138.7,114.1]
uppererror = [3.4,9.5,5.0,5.2,4.9,4.8,3.3,2.6]
lowererror = [3.4,9.5,5.0,5.2,4.9,4.8,3.3,2.6]
errs = [lowererror;uppererror]
res_tot = [scd1(ti) for ti in texp]
res = [ri[1] for ri in res_tot]
# Error computation    
wls_err_Cd1 = weighted_leastsquare(glcexpCd1,res,uppererror)
print("CD1 wls = ", wls_err_Cd1, "\n")
ls_err_Cd1 = leastsquare(glcexpCd1,res)
# Plotting    
plot!(scd1,vars=[(0,1)],label="Glucose(Cd1) simulations")
scatter!(texp,glcexpCd1,yerr = errs,marker = stroke(2, :green),
      xlabel="time (min)", ylabel= "Concentration: mg/dL",
      label="Cd1 exp results", legend=:top, markerstrokecolor = :auto)

    
# Parameters sets for the different hypothesis for the group CD2 at 26 Days
########################################################################    
# Control parameter set  ===> err 0.00403
#        p1    Gb     p2     p3    Ib     n       γ      h     p4      r
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 1.1: sensibility to interstitial insuline (r) ===> err :0.00304
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.05]
# hypothese 1.2: disponibility of insuline for X(t) (p3) ==> err: 0.00304
#  p0 = [0.01, 100.0, 0.79, 0.0354, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 2: degradation  insuline (n) ==> err : 0.00325
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.17, 0.0078, 65.0, 0.0170,  1.0]
# hypothese 3: interstitial insuline degradation (p2) ==> err :0.00312
#  p0 = [0.01, 100.0, 0.74, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothesis 4: γ (Insuline production acceleration) decreased  ==> err =  0.00641
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0170,  1.0]
# hypothesis 5: change p4 (stored insuline release rate)  ==> err = 0.00194      
  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 65.0, 0.0180,  1.0]
# hypothese 6: reponse Glucose threshold (h)  ==> err:  0.00402
#  p0 = [0.01, 100.0, 0.79, 0.0335, 06.0, 08.35, 0.0078, 66.0, 0.0170,  1.0]
########################################################################
# Simulation CD2 60 days    
prob_param = ODEProblem(f_param!,u0,tspan,p0)
scd2 = solve(prob_param, alg_hints = [:stiff])
# Experimental results    
glcexpCd2 = [97.9,141.0,186.1,196.0,196.6,172.5,132.5,111.4]
uppererror = [3.4,9.7,12.1,5.8,4.8,6.6,5.0,4.6]
lowererror = [3.4,9.7,12.1,5.8,4.8,6.6,5.0,4.6]
errs = [lowererror;uppererror]
res_tot = [scd2(ti) for ti in texp]
res = [ri[1] for ri in res_tot]
# Error computation    
wls_err_Cd2 = weighted_leastsquare(glcexpCd2,res,uppererror)
print("CD2 wls = ", wls_err_Cd2, "\n")
ls_err_Cd2 = leastsquare(glcexpCd2,res)
# Plotting    
plot!(scd2,vars=[(0,1)],label="Glucose (Cd2) simulations")
cur_colors = get_color_palette(:auto, plot_color(:white), 17)
scatter!(texp,glcexpCd2,yerr = errs,marker = stroke(2, cur_colors[5]),
       xlabel="time (min)", ylabel= "Concentration: mg/dL",
       label="Cd2 exp results", legend=:top, markerstrokecolor = :auto)

# Additional plot options    
plot!(size=(1280,960),legendfontsize=15, xtickfont=font(13), ytickfont=font(13), guidefont=font(18))
#  savefig("OGTT_60.pdf")
