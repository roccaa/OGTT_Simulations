# NEEDED MODULES -- This can take time at the first launch
# Code main be needed to run twice to plot correctly (shift-enter)
using DifferentialEquations
using Plots
using LaTeXStrings
# Load the backend (engine) for plotting based on javascript
plotlyjs()

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
    dose = 3729 # parameter K in paper
    return dose*1/(σ*sqrt(2*pi))*exp((-(t-μ)^2/(2*σ^2)))
end

# Parametric ODEs of te Glucose response model
# t : time
# u : variables
# p : parameters
# Here I call the guassian model of glucose input GlcGauss(t)
# And I saturate the variables to the positive cone
function f_param!(du,u,p,t)
 du[1] = -p[1]*(sat(u[1])-p[2])- p[10]*sat(u[2])*sat(u[1])+GlcGauss(t)
 du[2] = -p[3]*sat(u[2])+p[4]*(sat(u[3])-p[5])
 du[3] = -p[6]*sat(u[3])+p[7]*(sat(u[1])-p[8])*t + p[9]*u[4]
 du[4] = -p[9]*u[4]           
end   
  

# Initial conditions at 21j
u0 = [110.0,0.0,16.0,5950]
    
tspan = (0.0,120.0)

# #  A Parameter set from the MINMOD model as a reference
# #       p1  Gb    p2    p3     Ib  n   γ    h     |added parameters => |     p4    r
# p0 = [0.01,100.0,0.01,0.00001,8.0,0.1,0.01,100.0,                          0.07, 1.0]

  
# Fitted parameter set to the control group at 21 Days ===> err: 0.00311
#         p1    Gb    p2     p3     Ib     n      γ      h     p4     r
  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0] # err = 0.00311
###################################################
## Simulation of the Control experiment at 21Days
prob_param = ODEProblem(f_param!,u0,tspan,p0)  
sctrl = solve(prob_param, alg_hints = [:stiff])

## Experimental results    
texp = [0,10,20,30,45,60,90,120]
glcexp = [110.1,178.2,209.3,203.6,187.4,161.9,131.4,114.5]
uppererror = [1.4,4.0,3.3,2.6,2.8,2.4,1.9,1.5]
lowererror = [1.4,4.0,3.3,2.6,2.8,2.4,1.9,1.5]
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
    label="Control exp results", legend=:top, markerstrokecolor = :auto)

    
    
# Parameters sets for the different hypothesis for the group CD1 at 21 Days
########################################################################    
# Control parameter set  ===>   err = 0.0188
#         p1    Gb    p2     p3     Ib     n      γ      h     p4     r
  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0] 
# hypothese 1.1: sensibility to interstitial insuline (r) ===> err : 0.0144
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 0.88] 
# hypothese 1.2: disponibility of insuline for X(t) (p3) ==> err: 0.0144
#  p0 = [0.01, 100.0, 0.56, 0.0136, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0]
# hypothese 2: insuline degradation rate (n) ==> err : 0.0145
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 11.10, 0.0310, 85.0, 0.033, 1.0]
# hypothese 3: interstitial insuline degradation (p2) ==> err : 0.0140
# p0 = [0.01, 100.0, 0.65, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0]
# hypothesis 4: Slow phase 'acceleration' γ decreased  ==> err = 0.00566
 p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0258, 85.0, 0.033, 1.0] 
# hypothesis 5:  change p4 (stored insuline release rate)  ==> err = 0.0188
#   p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0]
# hypothese 6: reponse Glucose threshold (h)  ==> err: 0.0102
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 100.0, 0.033, 1.0] 
##############################################################################    
# Simulation CD1 21 Days    
prob_param = ODEProblem(f_param!,u0,tspan,p0) 
scd1 = solve(prob_param, alg_hints = [:stiff])
# Experimental results   
glcexpCd1 = [117.5,174.4,222.0,220.6,204.2,177.1,135.1,119.7]
uppererror = [1.8,2.7,3.8,4.1,3.7,2.9,2.2,1.5]
lowererror = [1.8,2.7,3.8,4.1,3.7,2.9,2.2,1.5]
errs = [lowererror;uppererror]
res_tot = [scd1(ti) for ti in texp]
res = [ri[1] for ri in res_tot]
# Error computation    
wls_err_Cd1 = weighted_leastsquare(glcexpCd1,res,uppererror)
print("CD1 wls = ", wls_err_Cd1, "\n")
ls_err_Cd1 = leastsquare(glcexpCd1,res)
# Plotting    
plot!(scd1,vars=[(0,1)],label="Glucose(Cd1) hyp 3 simu")
scatter!(texp,glcexpCd1,yerr = errs,marker = stroke(2, :green),
      xlabel="time (min)", ylabel= "Concentration: mg/dL",
      label="Cd1 exp results", legend=:top, markerstrokecolor = :auto)
    
# Parameters sets for the different hypothesis for the group CD2 at 21 Days
########################################################################    
# Control parameter set  ===> 0.0144 
#         p1    Gb    p2     p3     Ib     n      γ      h     p4     r
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0] 
# hypothese 1.1: sensibility to interstitial insuline (r) ==> err: 0.0082
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 0.70] 
# hypothese 1.2: disponibility of insuline for X(t) (p3) ==> err : 0.00821
#  p0 = [0.01, 100.0, 0.56, 0.0110, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0] 
# hypothese 2: degradation  insuline (n) ==> err : 0.00852
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 12.05, 0.0310, 85.0, 0.033, 1.0] 
# hypothese 3: interstitial insuline degradation (p2) ==> err: 0.00789
#  p0 = [0.01, 100.0, 0.81, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.033, 1.0] 
# hypothesis 4: Slow phase 'acceleration' γ decreased ==> err =  0.00433
 p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0215, 85.0, 0.033, 1.0]  
# hypothesis 5:  change p4 (stored insuline release rate)  ==> err = 0.01370
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 85.0, 0.029, 1.0]     
# hypothese 6: reponse Glucose threshold (h)  ==> err: 0.00740
#  p0 = [0.01, 100.0, 0.56, 0.0155, 10.0, 10.53, 0.0310, 115.0, 0.033, 1.0] 
############################################################################## 
# Simulation CD2 21 Days   
prob_param = ODEProblem(f_param!,u0,tspan,p0) 
scd2 = solve(prob_param, alg_hints = [:stiff])
# Experimental results    
glcexpCd2 = [105.0,196.9,226.8,243.7,220.7,191.7,151.6,119.2]
uppererror = [2.9,10.2,10.1,9.6,7.7,5.4,5.0,3.9]
lowererror = [2.9,10.2,10.1,9.6,7.7,5.4,5.0,3.9]
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

# Additional Plot options
plot!(size=(1280,960),legendfontsize=15, xtickfont=font(13), ytickfont=font(13), guidefont=font(18))