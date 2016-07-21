library(HTDoseResponseCurve)
context("Normalize synergy data")

test_that("synergy works", {
   
    pkg = "HTDoseResponseCurve"
    fn_data = system.file("extdata", "sample_synergy_ray.txt", package = pkg)
    syn=read.table( fn_data, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    head(syn)
    ds=create_synergy_dataset( syn$sample_type, syn$treatment, syn$treatment_2, 
         syn$concentration, syn$concentration_2, syn$value, 
         negative_control="DMSO")
    expect_equal( is_dataset_valid(ds), 2 )
    expect_equal( dim(syn)[1], length(ds$sample_type) )
    idx_l1_d = which( syn$sample_type=="line_1" & syn$treatment=="DMSO")
    mu_l1_d = mean(syn$value[idx_l1_d])
    
    idx_l1_6=which( syn$sample_type=="line_1" & syn$treatment=="drug_1" & 
               syn$concentration==6.25 & syn$concentration_2==0)
    expect_equal( round( syn$value[idx_l1_6] / mu_l1_d, 2 ), 
                  round( ds$value_normalized[idx_l1_6], 2))
    
    dose_SCH=c(0.1, 0.5, 1, 2, 4)
    eff_SCH = c(0.6701, 0.6289, 0.5577, 0.4550, 0.3755)
    dose_4HPR = c(0.1, 0.5, 1, 2)
    eff_4HPR = c(0.7666, 0.5833, 0.5706, 0.4934)
    eff_comb = c(0.6539, 0.4919, 0.3551, 0.2341)
    syn = data.frame( 
        treatment_1 = rep("SCH66336", 13),
        conc_1 = c( dose_SCH, rep(0, 4), dose_SCH[1:4]),
        treatment_2 = rep("4-HPR", 13),
        conc_2 = c( rep(0, 5), dose_4HPR, dose_4HPR ),
        values = c(eff_SCH, eff_4HPR, eff_comb ) )
   expect_error( create_synergy_dataset(
       sample_types = rep("sample_1", 13), 
       treatments_1 = syn$treatment_1,
       treatments_2 = syn$treatment_2,
       concentrations_1 = syn$conc_1,
       concentrations_2 = syn$conc_2,
       values = syn$values, negative_control = NA), 
       "treatment_1 and treatments_2 must be strings, not factors" )
   
} )