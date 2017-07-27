# SteadyCom python implementation
import scipy

def add_steadycom_variables_constraints(model,species_prefix_to_biomass_id,initial_guess=1.0):
    '''

    '''
    community_growth_rate = initial_guess
    biomass_vars = []
    for species_prefix in species_prefix_to_biomass_id.keys():
        # Create a biomass variable for each species and constrain to be positive
        species_biomass = model.problem.Variable(species_prefix+'_biomass')
        species_biomass_lb = model.problem.Constraint(
            species_biomass,name=species_prefix+'_biomass_bounds',lb=0,ub=1000)
        # hold a reference to the biomass variable for each species so that we can sum
        # it later in the total biomass constraint.
        biomass_vars.append(species_biomass)
        # couple flux through biomass with the community growth rate
        # NOTE: We assume that all fluxes are now abundance-scaled, meaning the flux through
        # biomass in the model is already abundance-scaled.
        community_coupling = model.problem.Constraint(
            model.reactions.get_by_id(species_prefix_to_biomass_id[species_prefix]).flux_expression - species_biomass,name=species_prefix+'_community_coupling', lb=0,ub=0)
        # add the constraints and variables generated so far to the model
        model.add_cons_vars([species_biomass,species_biomass_lb,community_coupling])
        # incorporate multiplication by the community growth rate by modifying the linear constraint on species_biomass
        # This needs to be done after the constraint is added to the model.
        model.constraints[species_prefix+'_community_coupling'].set_linear_coefficients({species_biomass:-1.0*community_growth_rate})
        # convert all reaction bounds to be abundance-scaled
        conversions_to_biomass = {}
        constraint_list = []
        for reaction in model.reactions:
            if reaction.id.find(species_prefix) >= 0:
                # for speed-up, add an empty constraint for each reaction, then set the linear coefficients directly
                # after constraints are added.
                constraint_lower = model.problem.Constraint(reaction.flux_expression - species_biomass, lb=0,name=reaction.id+'_lb_abundance_coupling')
                constraint_upper = model.problem.Constraint(reaction.flux_expression - species_biomass, ub=0,name=reaction.id+'_ub_abundance_coupling')
                model.add_cons_vars([constraint_lower,constraint_upper])
                model.constraints[reaction.id+'_lb_abundance_coupling'].set_linear_coefficients({reaction.forward_variable: 1,
                                                    reaction.reverse_variable: -1,
                                                    species_biomass: -1.0*reaction.lower_bound})

                model.constraints[reaction.id+'_ub_abundance_coupling'].set_linear_coefficients({reaction.forward_variable: 1,
                                                    reaction.reverse_variable: -1,
                                                    species_biomass: -1.0*reaction.upper_bound})
    # make individual biomass across all species sum to the total biomass
    total_biomass = model.problem.Variable('total_biomass')
    biomass_relations = model.problem.Constraint(
        total_biomass - sum(biomass_vars), lb=0, ub=0)
    # keeping the total biomass variable around can be useful for debugging/interpretation,
    # but the MaxBM problem doesn't constrain the total biomass, X0, to take a particular value (e.g. 1)
    model.add_cons_vars([total_biomass,biomass_relations])
    # update the objective to be total biomass, which is equal to the sum of species biomasses
    model.objective = total_biomass

def update_mu(model,species_prefix_to_biomass_id,new_guess):
    '''

    '''
    # For each species, update the community coupling constraint using the new_guess for mu.
    for species_prefix in species_prefix_to_biomass_id.keys():
        community_coupling = [constraint for constraint in model.constraints if constraint.name == species_prefix+'_community_coupling'][0]
        biomass_var = [model_variable for model_variable in model.variables if model_variable.name == species_prefix+'_biomass'][0]
        community_coupling.set_linear_coefficients({biomass_var:-1.0*new_guess})

def calculate_new_mu(old_mu,mu_lb,mu_ub,maxBM_solution,target_abundance=1.0,percent_change=0.001):
    if maxBM_solution >= target_abundance:
        mu_lb = old_mu
        change_coefficient = max([maxBM_solution/target_abundance, 1.0+percent_change])
        new_mu = change_coefficient*old_mu
    else:
        mu_ub = old_mu
        change_coefficient = max([maxBM_solution/target_abundance, 1.0-percent_change])
        new_mu = change_coefficient*old_mu
    return {'new_mu':new_mu,'new_lb':mu_lb,'new_ub':mu_ub}

def find_mu(model,mu_lb,mu_ub,target_abundance,species_prefix_to_biomass_id):
    '''
    Given a steadycom model and a range of mu values, find the species abundances,
    Xk, that satisfy f(mu) = X0.
    '''
    def abundances_minus_target_abundance(mu):
        #set mu in the model
        update_mu(model,species_prefix_to_biomass_id,mu)
        # calculate maxBM for the given community growth rate
        maxBM = model.optimize().f
        # we need to find the community growth rate within the interval
        # that satisfies maxBM = X0 (equivalent to X0 - maxBM)
        root_to_find = target_abundance - maxBM
        return root_to_find

    output = scipy.optimize.bisect(abundances_minus_target_abundance,\
                        mu_lb, mu_ub,)
    return output
