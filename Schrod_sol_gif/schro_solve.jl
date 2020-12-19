#Schrodinger Solver
using Plots

#FDTD - Finite-Difference Time-Domain

#Solver
function FDTD_solve(grid_step, x_lim, t_step, mass, pot_V, init_psi, n_step, bound_val=nothing)
    #Initialize
    grid_num = Int((x_lim[2]-x_lim[1])/grid_step) + 1
    x_space = LinRange(x_lim[1], x_lim[2], grid_num)
    psi_init_array = [init_psi(x) for x = x_space]
    pot_V_array = [pot_V(x) for x = x_space]

    psi_out = zeros(Complex{Float64}, (n_step, grid_num))
    psi_out[1, :] = psi_init_array
    function R2I_calc(x::Int, y::Int, grid_num::Int)
        if ((x == 1 && y == 1) || (x == grid_num && y == grid_num))
            return -1
        elseif (x == y)
            return -2
        elseif (abs(x-y) == 1)
            return 1
        else
            return 0
        end
    end
    diffM = [R2I_calc(x, y, grid_num) for x=1:grid_num, y=1:grid_num]
    psi_imag = imag.(psi_init_array)
    psi_real = real.(psi_init_array)

    planck = 1.054571817e-34
    c1 = (planck*t_step)/(2*mass*grid_step^2)
    c2 = t_step/planck
    for n in 2:n_step
        #Implementing update equations
        if (bound_val !== nothing)
            psi_real[1] = bound_val[1]
            psi_imag[grid_num] = bound_val[2]
        end
        psi_imag += c1*diffM*psi_real - c2*pot_V_array.*psi_real
        psi_real += -c1*diffM*psi_imag + c2*pot_V_array.*psi_imag
        psi_out[n, :] = psi_real.+psi_imag.*(1im)
    end
    return (x_space, pot_V_array, psi_out)
end




#Parameters and check for stability
n_steps = 50000
plank = 1.054571817e-34
mass = 2e-50
grid_step = 0.01
t_step = 1e-20
x_lim = [-5.0 20.0]
eV = 500
joules = 1.60217e-19*eV
pot_eV = 600
pot_max = 1.60217e-19*pot_eV
k = (2*(mass)*(joules)/(plank)^2)^(1/2)
stable = ((plank/(plank^2/(mass*grid_step^2) + pot_max/2)) > t_step, t_step-(plank/(plank^2/(mass*grid_step^2) + pot_max/2)))
not_alias = (10*mass < (plank^2)/(2*(grid_step)^2*joules), 10*mass-(plank^2)/(2*(grid_step)^2*joules))

#Potentials
function pot_well(x, bounds, height)
    if ((x <= bounds[1]) || (x >= bounds[2]))
        return height
    else
        return 0.0
    end
end
function pot_barr(x, bounds, height)
    if ((x >= bounds[1]) && (x <= bounds[2]))
        return height
    else
        return 0.0
    end
end
wave_fun(k, x) = cos(k*x) + sin(k*x)*(1im)
gaussian(x, mu, sig) = exp((-(1/2)*(x-mu)^2/(sig^2)))/(sig*(2*pi)^(1/2))
init_psi(x) = exp(-x^2+(1im)*17*x)*(2/pi)^(1/4)
pot_V_well(x) = pot_well(x, [-1.0, 1.0], pot_max)
pot_V_emp(x) = 0
pot_V_barr(x) = pot_barr(x, [5.0 5.25], pot_max)
pot_V_harm(x) = pot_max*(x/max(x_lim[1], x_lim[2]))^2
function pot_V_H(x)
    if (x == 0.0)
        return -pot_max
    else
        return (-(grid_step*pot_max)/abs(x) + pot_well(x, [-2.0, 2.0], pot_max))
    end
end

#Solution
#free_sol = FDTD_solve(grid_step, x_lim, t_step, mass, pot_V_emp, init_psi, n_steps)
#free_x_sol = free_sol[1]
#free_pot_sol = free_sol[2]
#free_real_sol = real.(free_sol[3])
#free_imag_sol = imag.(free_sol[3])
#free_prob_dense = free_real_sol.^2 + free_imag_sol.^2

#tunn_sol = FDTD_solve(grid_step, x_lim, t_step, mass, pot_V_barr, init_psi, n_steps)
#tunn_x_sol = tunn_sol[1]
#tunn_pot_sol = tunn_sol[2]./(1.60217e-19*pot_eV)
#tunn_real_sol = real.(tunn_sol[3])
#tunn_imag_sol = imag.(tunn_sol[3])
#tunn_prob_dense = tunn_real_sol.^2 + tunn_imag_sol.^2

#well_sol = FDTD_solve(grid_step, x_lim, t_step, mass, pot_V_well, init_psi, n_steps)
#well_x_sol = well_sol[1]
#well_pot_sol = well_sol[2]./(1.60217e-19*pot_eV)
#well_real_sol = real.(well_sol[3])
#well_imag_sol = imag.(well_sol[3])
#well_prob_dense = well_real_sol.^2 + well_imag_sol.^2

harm_sol = FDTD_solve(grid_step, x_lim, t_step, mass, pot_V_harm, init_psi, n_steps)
harm_x_sol = harm_sol[1]
harm_pot_sol = harm_sol[2]./(1.60217e-19*pot_eV)
harm_real_sol = real.(harm_sol[3])
harm_imag_sol = imag.(harm_sol[3])
harm_prob_dense = harm_real_sol.^2 + harm_imag_sol.^2

H_sol = FDTD_solve(grid_step, x_lim, t_step, mass, pot_V_H, init_psi, n_steps)
H_x_sol = H_sol[1]
H_pot_sol = H_sol[2]./(1.60217e-19*pot_eV)
H_real_sol = real.(H_sol[3])
H_imag_sol = imag.(H_sol[3])
H_prob_dense = H_real_sol.^2 + H_imag_sol.^2

#Plotting
#GIF-free
n_sec = 5
frame_rate = 30
df = Int(n_steps/(n_sec*frame_rate))
#free_anim = @animate for i = 1:df:n_steps
#    plot(free_x_sol, [free_real_sol[i,:] free_imag_sol[i,:] free_prob_dense[i,:] free_pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"], title = "Free Particle Solution", ylims = (-1, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
#well_anim = @animate for i = 1:df:n_steps
#    plot(well_x_sol, [well_real_sol[i,:] well_imag_sol[i,:] well_prob_dense[i,:] well_pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"], title = "Potential Well Solution", ylims = (-1, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
#tunn_anim = @animate for i = 1:df:n_steps
#    plot(tunn_x_sol, [tunn_real_sol[i,:] tunn_imag_sol[i,:] tunn_prob_dense[i,:] tunn_pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"], title = "Tunneling Solution", ylims = (-1, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
harm_anim = @animate for i = 1:df:n_steps
    plot(harm_x_sol, [harm_real_sol[i,:] harm_imag_sol[i,:] harm_prob_dense[i,:] harm_pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"], title = "QHO Solution", ylims = (-1, 1), ylabel = "Probability")
    println(i/n_steps*100)
end
H_anim = @animate for i = 1:df:n_steps
    plot(H_x_sol, [H_real_sol[i,:] H_imag_sol[i,:] H_prob_dense[i,:] H_pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"], title = "Hydrogen Atom Solution", ylims = (-1, 1), ylabel = "Probability")
    println(i/n_steps*100)
end

#free_anim_less = @animate for i = 1:df:n_steps
#    plot(free_x_sol, [free_prob_dense[i,:] free_pot_sol], label = ["ProbDense" "Potential"], title = "Free Particle Solution", ylims = (0, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
#well_anim_less = @animate for i = 1:df:n_steps
#    plot(well_x_sol, [well_prob_dense[i,:] well_pot_sol], label = ["ProbDense" "Potential"], title = "Potential Well Solution", ylims = (0, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
#tunn_anim_less = @animate for i = 1:df:n_steps
#    plot(tunn_x_sol, [tunn_prob_dense[i,:] tunn_pot_sol], label = ["ProbDense" "Potential"], title = "Tunneling Solution", ylims = (0, 1), ylabel = "Probability")
#    println(i/n_steps*100)
#end
harm_anim_less = @animate for i = 1:df:n_steps
    plot(harm_x_sol, [harm_prob_dense[i,:] harm_pot_sol], label = ["ProbDense" "Potential"], title = "QHO Solution", ylims = (0, 1), ylabel = "Probability")
    println(i/n_steps*100)
end
H_anim_less = @animate for i = 1:df:n_steps
    plot(H_x_sol, [H_prob_dense[i,:] H_pot_sol], label = ["ProbDense" "Potential"], title = "Hydrogen Atom Solution", ylims = (0, 1), ylabel = "Probability")
    println(i/n_steps*100)
end


#gif(free_anim, "free_schrod_sol.gif", fps = frame_rate)
#gif(well_anim, "well_schrod_sol.gif", fps = frame_rate)
#gif(tunn_anim, "tunn_schrod_sol.gif", fps = frame_rate)
gif(harm_anim, "harm_schrod_sol.gif", fps = frame_rate)
gif(H_anim, "h_schrod_sol.gif", fps = frame_rate)
#gif(free_anim_less, "free_schrod_sol_less.gif", fps = frame_rate)
#gif(well_anim_less, "well_schrod_sol_less.gif", fps = frame_rate)
#gif(tunn_anim_less, "tunn_schrod_sol_less.gif", fps = frame_rate)
gif(harm_anim_less, "harm_schrod_sol_less.gif", fps = frame_rate)
gif(H_anim_less, "h_schrod_sol_less.gif", fps = frame_rate)

#IMAGE
#plot(x_sol, [real_sol[i,:] imag_sol[i,:] prob_dense[i,:] pot_sol], label = ["Real" "Imag" "ProbDense" "Potential"])
#savefig("/home/6338_proj/plots/sol.png")
