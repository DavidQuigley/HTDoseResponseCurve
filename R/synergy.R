
#' Make two treatments consistently assigned to either treatment or treatment_2
#' 
#' Will raise an error if there are more than two unique treatments in the 
#' dataset. Applied across all hours and plate IDs.
#' 
#' @return A dataset with the same values as the dataset passed in as a 
#' parameter, but the treatment column will always be vehicle for treatment 1 or 
#' the value passed for treatment_1, and the treatment_2 column will always be 
#' the vehicle for treatment 2 or the value passed for treatment_2.
#' @param D dataset
#' @param treatment_1 String corresponding to the treatment to ensure is in 
#' the "treatment" column 
#' @param treatment_2 String corresponding to the treatment to ensure is in 
#' the "treatment_2" column 
#' @export
standardize_treatment_assigments = function(D, treatment_1, treatment_2){
    # confine treatment t1 to treatment
    # confine treatment t2 to treatment_2
    ds_type = is_dataset_valid(D, treatments=c(treatment_1, treatment_2))
    if( ds_type != 2 ){
        stop("This function can only be called on synergy datasets")   
    }
    treatments = unique( c( D$treatment[!D$is_negative_control], 
                            D$treatment_2[!D$is_negative_control_2] ) )
    if( length(treatments) != 2 ){
        stop(paste("There must be exactly two treatments in the dataset;",
                   "found:", paste(treatments, sep=",", collapse=",")))
    }
    for(i in 1:dim(D)[1] ){
        if( D$treatment[i]==treatment_2 | D$treatment_2[i]==treatment_1 ){
            tmp_conc = D$concentration[i]
            tmp_treat = D$treatment[i]
            tmp_is_neg = D$is_negative_control[i]
            tmp_neg = D$negative_control[i]
            D$concentration[i] = D$concentration_2[i]
            D$treatment[i] = D$treatment_2[i]
            D$is_negative_control[i] = D$is_negative_control_2[i]
            D$negative_control[i] = D$negative_control_2[i]
            D$concentration_2[i] = tmp_conc
            D$treatment_2[i] = tmp_treat
            D$is_negative_control_2[i] = tmp_is_neg
            D$negative_control_2[i] = tmp_neg
        }
    }
    D
}


split_synergy_datasets_for_CI = function( ds, treatment_1, treatment_2, 
                                          sample_type, hour ){
    is_treatment_1_vehicle = ds$is_negative_control 
    is_treatment_2_vehicle = ds$is_negative_control_2
    if( sum( is_treatment_1_vehicle )== 0 ){
        is_treatment_1_vehicle = ds$concentration==0
    }
    if( sum( is_treatment_2_vehicle )== 0 ){
        is_treatment_2_vehicle = ds$concentration_2==0
    }
    ds_1 = ds[ds$sample_type==sample_type & 
                  ds$hours==hour & 
                  ds$treatment==treatment_1 &
                  is_treatment_2_vehicle, ]
    ds_2 = ds[ds$sample_type==sample_type & 
                  ds$hours==hour & 
                  ds$treatment_2==treatment_2 &
                  is_treatment_1_vehicle, ]
    ds_c = ds[ds$sample_type==sample_type & 
                  ds$hours==hour & 
                  ds$treatment==treatment_1 &
                  ds$treatment_2==treatment_2 &
                  !is_treatment_1_vehicle & 
                  !is_treatment_2_vehicle, ]
    
    ds_1$conc_final = ds_1$concentration
    ds_2$conc_final = ds_2$concentration_2
    ds_c$conc_final = ds_c$concentration + ds_c$concentration_2
    list( a=ds_1, b=ds_2, ab=ds_c)
}

chou_synergy_helper = function( D, fu, fct ){
    # m is the slope of the median effect plot
    # measures the sigmoidicity of the dose effect curve; m=1 is hyperbolic
    fa = 1-fu
    fu[fu<0] = 0.001
    fa[fa<0] = 0.001
    y_m = rep(0, length(fu) )
    not_zero = !fa==0 & !fu==0
    y_m[ not_zero ] = log( fa[not_zero] / fu[not_zero] )
    
    x_m = rep(0, length(D))
    not_zero = D != 0
    x_m[not_zero] = log( D[not_zero] )

    m = as.numeric(stats::coefficients( stats::lm( y_m ~ x_m ) )[2])
    
    # Dm is the IC50 for this curve
    data = data.frame(fu, D)
    ed = drc::ED( drc::drm(fu~D, data=data, fct=fct), c(50), display=FALSE )
    Dm = ed[1]
    Dm_stderr = ed[2]
    list( D=D, fa=fa, fu=fu, y_m=y_m, x_m=x_m, m=m, Dm=Dm, Dm_stderr=Dm_stderr )
}

#' Calculate statistics needed for Chou Synergy
#' Curve fitting is performed by the \code{drm()} function in the \code{drc} 
#' library. To fit the curve, you need to select a non-linear function. To 
#' estimate the slope, upper asymptote, lower asymptote, and EC50, pass 
#' drc::LL.4(). To fix the lower asymptote at 1 and estimate the other 
#' parameters, pass drc::LL.3(). To fix the upper asympotote at 1 and the lower 
#' asymptote at 0, pass dcr::LL.2. For a list of available functions, see 
#' \code{drc::getMeanFunctions()}. 
#'
#' Ting-Chao Chou Pharmacological Reviews 2006
#' 
#' @return A list with the Chou Synergy helper values for treatment_1, 
#' treatment_2, the combination of treatment_12, and the combination index 
#' data frame (itself containing columns for Fa, Fu, and CI).
#' @param ds dataset
#' @param sample_type sample type in ds
#' @param treatment_1 treatment in ds
#' @param treatment_2 treatment in ds
#' @param hour hour in ds
#' @param fct Non-linear function to fit, e.g. drc::LL.3(). See summary.
#' @param summary_method mean or median
#' @export
chou_synergy = function( ds, sample_type, treatment_1, treatment_2, hour, 
                         fct, summary_method ){
    # m is the slope of the median effect plot
    # measures the sigmoidicity of the dose effect curve; m=1 is hyperbolic
    
    ds_type = is_dataset_valid(ds, 
                               treatments = c(treatment_1, treatment_2),
                               hour=hour)
    
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    
    if( ds_type != 2 ){ stop(paste("dataset ds must be a synergy dataset with",
                              "columns called treatment_2 and concentration_2"))   
    }
    
    a_b_ab = split_synergy_datasets_for_CI( ds, treatment_1, treatment_2, 
                                            sample_type, hour )
    ds_1 = a_b_ab$a
    ds_2 = a_b_ab$b
    ds_c = a_b_ab$ab
    
    proportion_1= ds_c$concentration/(ds_c$concentration + ds_c$concentration_2)
    proportion_2 = 1-proportion_1
    get_mu = function(x){ Fa=mean(x$value_normalized, na.rm=TRUE) }
    get_median = function(x){ Fa=stats::median(x$value_normalized, na.rm=TRUE) }
    
    if( summary_method=="mean" ){
        Dsum1 = plyr::ddply(ds_1, c("conc_final"), get_mu )
        Dsum2 = plyr::ddply(ds_2, c("conc_final"), get_mu )
        Dsumc = plyr::ddply(ds_c, c("conc_final"), get_mu )
    }else{
        Dsum1 = plyr::ddply(ds_1, c("conc_final"), get_median )
        Dsum2 = plyr::ddply(ds_2, c("conc_final"), get_median )
        Dsumc = plyr::ddply(ds_c, c("conc_final"), get_median )
    }
    names(Dsum1) = c("D", "Fu")
    names(Dsum2) = c("D", "Fu")
    names(Dsumc) = c("D", "Fu")
    cs_1 = chou_synergy_helper( Dsum1$D, Dsum1$Fu, fct )
    cs_2 = chou_synergy_helper( Dsum2$D, Dsum2$Fu, fct )
    cs_c = chou_synergy_helper( Dsumc$D, Dsumc$Fu, fct )
    
    Fa = cs_c$fa
    Fu = 1-cs_c$fa
    
    CI_1 = (cs_c$D * proportion_1) / (cs_1$Dm * (Fa/Fu)^(1/cs_1$m) )
    CI_2 = (cs_c$D * proportion_2) / (cs_2$Dm * (Fa/Fu)^(1/cs_2$m) )
    CI = data.frame( Fa=Fa, Fu=Fu, CI=CI_1+CI_2 )
    L=list( treatment_1 = cs_1, 
            treatment_2 = cs_2, 
            treatment_12 = cs_c, 
            CI=CI)
}


#' Plot single-value Chou Combination Index at IC50
#' 
#' @param CS chou statistics calculated by \code{\link{chou_synergy}}
#' @param proportion_1 proportion of dose 1 vs. dose 2, numeric between 0 and 1
#' @return Combination Index
#' @export
chou_synergy_CI_median = function( CS, proportion_1 ){
    D1 = CS$treatment_12$Dm * proportion_1
    D2 = CS$treatment_12$Dm * (1-proportion_1)
    Dm1 = CS$treatment_1$Dm
    Dm2 = CS$treatment_2$Dm
    (D1/Dm1) + (D2/Dm2)
}

#' Construct confidence intervals for observed effects at combination doses 
#' having observed effects.
#' 
#' If the dataset does not have negative controls, measurements at 
#' concentration zero are assumed to be empty for a given treatment.
#' Adapted directly from code published in 
#' Lee and Kong Statistics in Biopharmaceutical Research 2012
#' 
#' @return A list with the interaction index, lower, and upper confidence 
#' intervals for fix effects points in E
#' @param ds dataset
#' @param sample_type sample type in ds
#' @param treatment_1 treatment in ds
#' @param treatment_2 treatment in ds
#' @param hour hour in ds. Default 0. 
#' @param alpha 1-alpha is the size of the confidence intervals, default 0.05
#' @param summary_method mean or median
#' @references Lee & Kong Statistics in Biopharmaceutical Research 2012
#' @export
synergy_interaction_CI = function(ds, sample_type, 
                                  treatment_1, treatment_2,
                                  hour=0, alpha=0.05,
                                  summary_method="mean"){
    
    ds_type = is_dataset_valid(ds, treatments=c(treatment_1, treatment_2),
                               hour=hour)
    
    if(! (summary_method=="mean" | summary_method=="median" ) ){
        stop("parameter summary_method must be either mean or median")   
    }
    if( summary_method=="mean" ){
        sumfunc = function(po){ data.frame( 
            value=mean(po$value_normalized, na.rm=TRUE)) }
    }
    else{
        sumfunc = function(po){ data.frame( 
            value=stats::median(po$value_normalized, na.rm=TRUE)) }
    }
    
    a_b_ab = split_synergy_datasets_for_CI( ds, treatment_1, treatment_2, 
                                            sample_type, hour )
    ds_1 = a_b_ab$a
    ds_2 = a_b_ab$b
    ds_c = a_b_ab$ab
    
    if( dim(ds_1)[1]==0 )
        stop("No samples meet the criteria for treatment_1")
    if( dim(ds_2)[1]==0 )
        stop("No samples meet the criteria for treatment_2")
    if( dim(ds_c)[1]==0 )
        stop("No samples meet the criteria for combined treatments")
    
    E_ci = seq(from=0.01, to=1, by=0.01)
    E_combined = plyr::ddply( ds_c, c("conc_final"), 
                              function(po){ data.frame( 
                                  mu=mean(po$value_normalized, na.rm=TRUE)) 
                              } )$mu
    if( sum(E_combined>1)>0 ){
        warning("One or more effects were greater than 1, setting to 0.999")   
    }
    E_combined[ E_combined >= 1 ] = 0.999
    obs_plus_ci = data.frame( is_obs = c( rep( TRUE, length(E_combined)), 
                                          rep( FALSE, length(E_ci)) ),
                              effect = c(E_combined, E_ci ) )
    obs_plus_ci=obs_plus_ci[order(obs_plus_ci$effect),]
    
    E = obs_plus_ci$effect
    
    ds_1$conc_final = ds_1$concentration
    ds_2$conc_final = ds_2$concentration_2
    ds_c$conc_final = ds_c$concentration + ds_c$concentration_2
    E1 = plyr::ddply( ds_1, c("conc_final"), sumfunc )
    e1 = E1$value
    d1 = E1$conc_final
    
    E2 = plyr::ddply( ds_2, c("conc_final"), sumfunc)
    e2 = E2$value
    d2 = E2$conc_final
    
    E12 = plyr::ddply( ds_c, c("conc_final"), sumfunc)
    e12 = E12$value
    d12 = E12$conc_final
    
    e1[e1>1] = 0.999
    e2[e2>1] = 0.999
    e12[e12>1] = 0.999
    
    d2.d1 = ds_c$concentration_2 / ds_c$concentration
    if( length(unique(d2.d1))>1 ){
        stop("not all combinations were made at the same ratio")
    }
    d2.d1 = d2.d1[1]
    
    lm1 =  stats::lm( log(e1/(1-e1)) ~ log(d1) )
    lm2 =  stats::lm( log(e2/(1-e2)) ~ log(d2) )
    lm12 = stats::lm( log(e12/(1-e12)) ~ log(d12) )
    dm1 =  exp( -summary(lm1)$coef[1,1] / summary(lm1)$coef[2,1] )
    dm2 =  exp( -summary(lm2)$coef[1,1] / summary(lm2)$coef[2,1] )
    dm12 = exp( -summary(lm12)$coef[1,1] / summary(lm12)$coef[2,1] ) 
    
    Dx1 = dm1*(E/(1-E))^(1/summary(lm1)$coef[2,1])
    Dx2 = dm2*(E/(1-E))^(1/summary(lm2)$coef[2,1])
    dx12 = dm12*(E/(1-E))^(1/summary(lm12)$coef[2,1])
    iix = ( dx12 / (1+d2.d1) ) / Dx1 + (dx12 * d2.d1 / (1+d2.d1) ) / Dx2
    lm1.s = summary(lm1)
    lm2.s = summary(lm2)
    lm12.s = summary(lm12)
    c1 = 1.0 / lm1.s$coef[2,1]^2 * lm1.s$coef[1,2]^2
    temp = - mean(log(d1)) * lm1.s$coef[2,2]^2
    c1 = c1 + 2.0*(log(E/(1-E)) - lm1.s$coef[1,1]) / lm1.s$coef[2,1]^3 * temp
    c1 = c1+(log(E/(1-E))-lm1.s$coef[1,1])^2/lm1.s$coef[2,1]^4*lm1.s$coef[2,2]^2
    
    c2 = 1.0 / lm2.s$coef[2,1]^2 * lm2.s$coef[1,2]^2
    temp = - mean(log(d2)) * lm2.s$coef[2,2]^2
    c2 = c2 + 2.0 * (log(E/(1-E)) - lm2.s$coef[1,1]) / lm2.s$coef[2,1]^3 * temp
    c2 = c2+(log(E/(1-E))-lm2.s$coef[1,1])^2/lm2.s$coef[2,1]^4*lm2.s$coef[2,2]^2
    
    c12 = 1.0 / lm12.s$coef[2,1]^2 * lm12.s$coef[1,2]^2
    temp = - mean( log(d12) ) * lm12.s$coef[2,2]^2
    c12 = c12 + 2.0 * (log(E/(1-E)) - lm12.s$coef[1,1]) /lm12.s$coef[2,1]^3*temp
    c12 = c12+(log(E/(1-E))-lm12.s$coef[1,1])^2 / 
              lm12.s$coef[2,1]^4*lm12.s$coef[2,2]^2
    
    var.ii =( (dx12 / Dx1)^2 * c1 + ( dx12 * d2.d1 / Dx2 )^2 * 
                  c2+( 1.0 / Dx1 + d2.d1 / Dx2 )^2 * dx12^2 * c12) / (1+d2.d1)^2 
    t975 = stats::qt( 1-alpha/2, length(d1) + length(d2) + length(d12) - 6 )
    iix.lo1 = iix * exp( -t975 * var.ii^0.5 / iix )
    iix.hi1 = iix * exp(  t975 * var.ii^0.5 / iix )
    
    iix.lo1[ is.nan(iix.lo1) ] = NA
    iix.hi1[ is.nan(iix.hi1) ] = NA
    list(interaction_index = iix, 
         cl_lower = iix.lo1, 
         cl_upper = iix.hi1,
         effects = obs_plus_ci$effect,
         is_obs = obs_plus_ci$is_obs)
}


#' Construct synergy matrix of values expected from Bliss independence those 
#' actually observed
#' 
#' @param D dataset
#' @param treatment_1 treatment in D
#' @param treatment_2 treatment in D
#' @examples 
#' samples = rep("s1", 16)
#' t1 = rep( c("DMSO", "d1", "d1", "d1"), 4)
#' t2 = c( rep( "DMSO", 4), rep("d2", 12) )
#' c1 = rep( c(0, 50, 100, 200), 4)
#' c2 = c(0,0,0,0, 50,50,50,50, 100,100,100,100, 200,200,200,200)
#' value_ind=c(1,0.8,0.7,0.6,0.8,0.7,0.6,0.5,0.7,0.6,0.5,0.4,0.6,0.5,0.4,0.3)
#' value_syn=c(1,0.8,0.7,0.6,0.8,0.8,0.5,0.2,0.7,0.2,0.1,0.05,0.6,0.1,0.05,0.01)
#' 
#'  DS_i=create_synergy_dataset( sample_types=samples, treatments_1=t1,
#'                          treatments_2=t2, concentrations=c1, 
#'                          concentrations_2=c2, values=value_ind, 
#'                          negative_control="DMSO")
#'  DS_s=create_synergy_dataset( sample_types=samples, treatments_1=t1,
#'                          treatments_2=t2, concentrations=c1, 
#'                          concentrations_2=c2, values=value_syn, 
#'                          negative_control="DMSO")
#'                          
#' b_ind = synergy_bliss(DS_i, "d1", "d2")
#' b_syn = synergy_bliss(DS_s, "d1", "d2")
#' layout(matrix(1:2,1,2))
#' plot_color_grid(b_ind$excess, color_bounds = c(-1,1))
#' plot_color_grid(b_syn$excess, color_bounds = c(-1,1))
#' @return A list with two matrixes, one for the values expected under Bliss 
#' independence, a + b - (a*b), and one for what was actually observed. 
#' 
#' @export
#' 
synergy_bliss = function( D, treatment_1, treatment_2 ){
    if( length(D$treatment==treatment_1)==0 ){
        stop(paste("value",treatment_1,"not found in column treatment"))
    }
    if( length(D$treatment_2==treatment_2)==0 ){
        stop(paste("value",treatment_2,"not found in column treatment_2"))
    }    
    vehicle_1 = unique(D$negative_control[D$treatment==treatment_1] )
    vehicle_2 = unique(D$negative_control[D$treatment_2==treatment_2] )
    D = D[D$treatment %in% c(treatment_1, vehicle_1) |
              D$treatment_2 %in% c(treatment_2, vehicle_2) ,]
    
    D$value_normalized = 1-D$value_normalized
    
    cols = sort(unique(D$concentration[D$treatment==treatment_1]))
    rows = sort(unique(D$concentration_2[D$treatment_2==treatment_2]))
    bliss_expected = matrix(NA, ncol=length(cols), nrow=length(rows))
    bliss_observed = matrix(NA, ncol=length(cols), nrow=length(rows))
    dimnames(bliss_expected)[[1]] = rows
    dimnames(bliss_expected)[[2]] = cols
    dimnames(bliss_observed)[[1]] = rows
    dimnames(bliss_observed)[[2]] = cols
    for(c_idx in 1:length(cols)){
        for(r_idx in 1:length(rows)){
            conc_c = cols[c_idx]
            conc_r = rows[r_idx]
            a = D$value_normalized[ D$treatment==treatment_1 & 
                                     D$concentration==conc_c & 
                                     D$treatment_2 == vehicle_2 ]
            b = D$value_normalized[ D$treatment_2==treatment_2 & 
                                         D$concentration_2==conc_r & 
                                         D$treatment == vehicle_1 ]
            ab = D$value_normalized[ D$treatment==treatment_1 & 
                                         D$treatment_2==treatment_2 & 
                                         D$concentration==conc_c & 
                                         D$concentration_2==conc_r] 
            a = mean(a, na.rm=TRUE)
            b = mean(b, na.rm=TRUE)
            ab = mean( ab, na.rm=TRUE)
            if(is.nan(a)){ a=0 }
            if(is.nan(b)){ b=0 }
            bliss_expected[r_idx, c_idx] = a + b - (a*b)
            bliss_observed[r_idx, c_idx] = ab
        }
    }
    list( expected = bliss_expected, 
          observed = bliss_observed,
          excess = round( bliss_observed-bliss_expected, 2 ) )
}