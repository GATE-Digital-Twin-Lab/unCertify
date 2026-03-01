#----------------------------2-D cases----------------------------
#--------------------Cruise case----------------------------------
## ttw_crs = sc.thrust_to_weight_cruise(
#     wing_ld_pa, cruisealt_m, cruisespeed_mps, cd_min_clean, aspect_ratio)

def case_1(sc, wing_ld_pa, cruisealt_m, cruisespeed_mps, vars, AA_vars, base_params):
    ttw = sc.thrust_to_weight_cruise(
        wing_ld_pa, 
        cruisealt_m,
        cruisespeed_mps,
        AA_vars[3],
        AA_vars[5]
    )

    return {
        "f_target": ttw,
        "vars": [vars[3], vars[5]],
        "V_base": [base_params[3], base_params[5]],
        "param_indices": [0, 3, 5] 
    }

#--------------------Take off case----------------------------------

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_2(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        AA_vars[1],
        mid_base_params[2],
        groundrun_m,
        mid_base_params[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[1]],
        "V_base": [base_params[0], base_params[1]],
        "param_indices": [1, 0, 1] 
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_3(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        mid_base_params[0],
        AA_vars[1],
        AA_vars[2],
        groundrun_m,
        mid_base_params[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[1], vars[2]],
        "V_base": [base_params[1], base_params[2]],
        "param_indices": [1, 1, 2] 
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_4(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        mid_base_params[0],
        mid_base_params[1],
        AA_vars[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[2], vars[4]],
        "V_base": [base_params[2], base_params[4]],
        "param_indices": [1, 2, 4] 
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_5(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        mid_base_params[1],
        AA_vars[2],
        groundrun_m,
        mid_base_params[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[2]],
        "V_base": [base_params[0], base_params[2]],
        "param_indices": [1, 0, 2] 
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_6(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        mid_base_params[1],
        mid_base_params[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[4]],
        "V_base": [base_params[0], base_params[4]],
        "param_indices": [1, 0, 4] 
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_7(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        mid_base_params[0],
        AA_vars[1],
        mid_base_params[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[1], vars[4]],
        "V_base": [base_params[1], base_params[4]],
        "param_indices": [1, 1, 4] 
    }


#----------------------------3-D cases----------------------------
# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_8(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        AA_vars[1],
        AA_vars[2],
        groundrun_m,
        mid_base_params[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[1], vars[2]],
        "V_base": [base_params[0], base_params[1], base_params[2]],
        "param_indices": [1, 0, 1, 2]
    }


# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_9(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        AA_vars[1],
        mid_base_params[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[1], vars[4]],
        "V_base": [base_params[0], base_params[1], base_params[4]],
        "param_indices": [1, 0, 1, 4]
    }


# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_10(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        mid_base_params[0],
        AA_vars[1],
        AA_vars[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[1], vars[2], vars[4]],
        "V_base": [base_params[1], base_params[2], base_params[4]],
        "param_indices": [1, 1, 2, 4]
    }

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_11(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        mid_base_params[1],
        AA_vars[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[2], vars[4]],
        "V_base": [base_params[0], base_params[2], base_params[4]],
        "param_indices": [1, 0, 2, 4]
    }


#----------------------------4-D cases----------------------------

# ttw_to = sc.thrust_to_weight_take_off(
#     wing_ld_pa, cl_max_to, cl_to, cd_to, groundrun_m, mu_r) 

def case_12(sc, wing_ld_pa, vars, AA_vars, base_params, mid_base_params, groundrun_m):
    ttw = sc.thrust_to_weight_take_off(
        wing_ld_pa,
        AA_vars[0],
        AA_vars[1],
        AA_vars[2],
        groundrun_m,
        AA_vars[4]
    )

    return {
        "f_target": ttw,
        "vars": [vars[0], vars[1], vars[2], vars[4]],
        "V_base": [base_params[0], base_params[1], base_params[2], base_params[4]],
        "param_indices": [1, 0, 1, 2, 4]
    }