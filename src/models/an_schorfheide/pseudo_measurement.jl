"""
```
pseudo_measurement{T<:AbstractFloat}(m::AnSchorfheide{T};
    apply_altpolicy::Bool = false)
```

Assign pseudo-measurement equation (a linear combination of states):

```
X_t = ZZ_pseudo*S_t + DD_pseudo
```
"""
function pseudo_measurement{T<:AbstractFloat}(m::AnSchorfheide{T}; apply_altpolicy::Bool = false)

    endo = m.endogenous_states
    pseudo_names = [:y_t, :π_t, :z_t, :NominalFFR, :RealFFR]

    # Map pseudoobservables to indices
    pseudo_inds = OrderedDict{Symbol,Int}()
    for (i,k) in enumerate(pseudo_names)
        pseudo_inds[k] = i
    end

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(length(pseudo), n_states_augmented(m))
    DD_pseudo = zeros(length(pseudo))

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo_inds[:y_t],endo[:y_t]] = 1.
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    ## π_t
    ZZ_pseudo[pseudo_inds[:π_t],endo[:π_t]] = 1.
    DD_pseudo[pseudo_inds[:π_t]]            = 100*(m[:π_star]-1);

    pseudo[:π_t].name = "Inflation"
    pseudo[:π_t].longname = "Inflation"
    pseudo[:π_t].rev_transform = quartertoannual

    ## z_t
    ZZ_pseudo[pseudo_inds[:z_t], endo[:z_t]] = 1.
    pseudo[:z_t].name     = "z_t"
    pseudo[:z_t].longname = "z_t"

    ## Nominal FFR
    ZZ_pseudo[pseudo_inds[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo_inds[:NominalFFR]] = m[:π_star] + m[:rA] + 4.0*m[:γ_Q]
    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"

    ## Real FFR
    ZZ_pseudo[pseudo_inds[:RealFFR], endo[:R_t]] = 1.
    ZZ_pseudo[pseudo_inds[:RealFFR], endo[:π_t]] = -4.
    DD_pseudo[pseudo_inds[:RealFFR]] = m[:π_star] + m[:rA] + 4.0*m[:γ_Q] - 4.0*(100. * (m[:π_star] - 1.))
    pseudo[:RealFFR].name     = "Real FFR"
    pseudo[:RealFFR].longname = "Real FFR at an annual rate"

    # Collect indices and transforms
    return pseudo, PseudoObservableMapping(pseudo_inds, ZZ_pseudo, DD_pseudo)
end