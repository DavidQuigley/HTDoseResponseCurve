library(HTDoseResponseCurve)
library(drc)
context("fit and summary methods")

test_that("fit works", {
   
    # # line 1 has measurements at 24 and 48 hours while line 2 has only 48 hours
    # ds1 = create_dataset( 
    #     sample_types= c("line1","line1","line1","line1","line1","line1",
    #                     "line2","line2","line2"),
    #     treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2",
    #                    "DMSO","drug1","drug2"),
    #     concentrations = c(0, 100, 200, 0, 100, 200, 
    #                        0, 100, 200),
    #     hours = c(24, 24, 24, 48, 48, 48, 
    #               48, 48, 48),
    #     values = c(100, 90, 20, 99, 80, 15,
    #                100, 89, 87), 
    #     plate_id = "plate_1",
    #     negative_control = "DMSO")
    # 
    # fits = fit_statistics(ds1, drc::LL.4() )
    # 
    # ds2 = create_dataset(
    #     sample_types= c("line1","line1","line1","line1",
    #                     "line1","line1","line1","line1",
    # 
    #                     "line2","line2","line2","line2",
    #                     "line2","line2","line2","line2",
    # 
    #                     "line3","line3","line3","line3",
    #                     "line3","line3","line3","line3"),
    #     treatments = rep(c("DMSO","drug1","drug1","drug1"), 6),
    #     concentrations = c(0, 50, 100, 200,
    #                        0, 50, 100, 200,
    # 
    #                        0, 50, 100, 200,
    #                        0, 50, 100, 200,
    # 
    #                        0, 50, 100, 200,
    #                        0, 50, 100, 200),
    #     hours = rep(0, 24),
    #     values = c(100, 90, 85, 72,
    #                99, 92, 83, 79,
    # 
    #                97, 80, 40, 15,
    #                95, 78, 38, 13,
    # 
    #                96, 60, 35, 10,
    #                94, 58, 33, 12),
    #     plate_id = "plate_1",
    #     negative_control = "DMSO")
    # 
    # fit=fit_DRC( ds2,
    #              sample_types = c("line1", "line2", "line3"),
    #              treatments = c("drug1"),
    #              hour=0, fct=drc::LL.3() )
    # fpc=HT_fit_plot_colors( fit )
    # plot( fit,
    #       xlim=c(0, 1e5), ylim=c(0, 1.2),
    #       lwd=3,
    #       main="test",
    #       ylab="surviving fraction",
    #       xlab="nM")
    # legend( 2, 0.5, legend = get.split.col(fpc$condition, "_|_", first = TRUE),
    #         fill=fpc$col, bty="n")
    # 
    # fits = fit_statistics(ds2, drc::LL.4() )
    # 
} )
