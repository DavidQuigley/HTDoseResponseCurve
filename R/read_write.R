#' Read plate maps from Incucyte plate map XML extract
#'
#' A Plate map describes the cell lines, treatment, and treatment concentrations
#' applied to each well. This file is a text file in the XML format.
#' 
#' @param path_to_file The complete path to a plate map file.
#' @param attribute_value The XML attribute to extract for treatment and sample 
#' names, one of {"description", "displayName"}. Defaults to "description".
#' @param max_treatments_per_well A number indicating the maximum number of 
#' treatments to expect in a single well. For single-treatment experiments, pass 
#' 1. For synergy experiments where two treatments may be combined in a single 
#' well, pass 2. Defaults to 1.
#' @return A list of matrixes with the same dimension as the plate map. The 
#' matrixes are concentration, treatment, sample_type, seeding_density, and 
#' passage
#' @import XML
#' @examples 
#' pkg = "HTDoseResponseCurve"
#' fn_map = system.file("extdata", "sample_data_384_platemap.txt",package=pkg)
#' plate_map = read_platemap_from_Incucyte_XML( fn_map )
#' @export
read_platemap_from_Incucyte_XML = function( path_to_file, 
                                            attribute_value="description",
                                            max_treatments_per_well = 1){
    if( !file.exists(path_to_file) ){
        stop(paste("Cannot open a file at",path_to_file))   
    }
    if( !is.numeric(max_treatments_per_well)){
        stop("parameter max_treatments_per_well must be either 1 or 2")   
    }
    if( max_treatments_per_well!=1 & max_treatments_per_well!=2){
        stop( "parameter max_treatments_per_well must be either 1 or 2")   
    }
    if( attribute_value != "description" & attribute_value != "displayName"){
        stop(paste("attribute_value parameter must be either description or",
                   "displayName"))
    }
    data = XML::xmlParse(path_to_file)
    x = XML::xmlToList(data)
    n_wells = length( x$wellStore$wells )
    pm = create_empty_plate_map( n_wells,
                                max_treatments_per_well=max_treatments_per_well)
    for(i in 1 : n_wells ){
        if( ".attrs" %in% names(x$wellStore$wells[i]$well)){
            well_attr = x$wellStore$wells[i]$well$.attrs
            rr = as.numeric( well_attr[ names(well_attr)=="row" ] ) + 1
            cc = as.numeric( well_attr[ names(well_attr)=="col" ] ) + 1
            well_items = x$wellStore$wells[i]$well$items
            compounds_seen = 0
            for( j in 1: length( well_items ) ){
                well_ref = well_items[j]$wellItem$referenceItem
                reference_type = as.character(well_ref[names(well_ref)=="type"])
                if( reference_type=="CellType" ){
                    idx_displayname=which(names(well_ref)==attribute_value)
                    pm$sample_type[rr,cc]=as.character(well_ref[idx_displayname])
                    well_attrs = well_items[j]$wellItem$.attrs
                    idx_density = which(names(well_attrs)=="seedingDensity")
                    idx_passage = which(names(well_attrs)=="passage")
                    pm$density[rr,cc] = as.numeric( well_attrs[ idx_density ] )
                    pm$passage[rr,cc] = as.numeric( well_attrs[ idx_passage ] )
                }else if( reference_type=="Compound" ){
                    idx_display = which(names(well_ref)==attribute_value)
                    if( compounds_seen==0 ){
                        pm$treatment[rr,cc]=as.character(well_ref[idx_display])
                        well_attrs = well_items[j]$wellItem$.attrs
                        idx_conc = which(names(well_attrs)=="concentration")
                        pm$concentration[rr,cc]=as.numeric(well_attrs[idx_conc])
                        compounds_seen = 1
                    }else{
                        if( max_treatments_per_well == 1 ){
                            warning(paste("parameter max_treatments_per_well",
                             "was 1, but detected second treatment in at least",
                             "one well. Is this a synergy experiment?"))
                        }
                        pm$treatment_2[rr,cc]=as.character(well_ref[idx_display])
                        well_attrs = well_items[j]$wellItem$.attrs
                        idx_conc = which(names(well_attrs)=="concentration")
                        pm$concentration_2[rr,cc]=as.numeric(well_attrs[idx_conc])
                        compounds_seen = 2
                    }
                    
                }
            }
        }
    }
    pm$treatment[pm$treatment==""] = NA
    pm$sample_type[pm$sample_type==""] = NA
    pm$density[pm$density==""] = NA
    pm$passage[pm$passage==""] = NA
    if( "concentration_2" %in% names(pm) ){
        pm$concentration_2[pm$concentration_2==""] = NA   
    }
    if( "treatment_2" %in% names(pm) ){
        pm$treatment_2[pm$treatment_2==""] = NA   
    }
    pm
} 


#' Read plate maps into a list of three data frames
#'
#' A Plate map describes the cell lines, treatment, and treatment concentrations
#' applied to each well. As an alternate to exporting a plate map into 
#' XML directly from the instrument, \code{read_platemap_from_excel} loads plate
#' map data from an Excel file. The input file must contain cells in column A
#' indicating the values for treatment concentration in uM (e.g. 0.2 for 200 
#' nM), treatment (e.g. the drug in a well), and sample type (e.g. the 
#' identifier for a particular cell line in a well). These cells should be 
#' one cell up and to the left of plate maps themselves.
#' 
#' @param filename Excel file generated by the user
#' @param sheet_num Index of the Excel sheet to read from filename. Defaults to 
#' 1.
#' @param number_of_wells The number of wells in each plate; must be one of 
#' 6, 12, 24, 96, 384.
#' @param concentration_identifier Text in a cell in column 1 that signals the 
#' concentration map will follow. Defaults to "concentration".
#' @param treatment_identifier text in a cell in column 1 that signals the 
#' treatment_identifier map will follow. Defaults to "treatment".
#' @param sample_identifier text in a cell in column 1 that signals the 
#' sample_identifier map will follow. Defaults to "cell line".
#' @param na.strings text in a cell that is interpreted as missing data. 
#' Defaults to the strings "NA" or an empty cell.
#' @return A data frame where columns are data, timestamp, plate_id, hour.
#' @import readxl
#' @examples
#' pkg = "HTDoseResponseCurve"
#' fn_map = system.file("extdata", "sample_data_96_platemap.xlsx",package=pkg)
#' plate_map = read_platemap_from_excel( fn_map, number_of_wells=96 )
#' @export
read_platemap_from_excel = function( filename, sheet_num=1, number_of_wells,
                                     concentration_identifier="concentration",
                                     treatment_identifier="treatment",
                                     sample_identifier="cell line",
                                     na.strings=c("", "NA") ){
    
    xl = data.frame( readxl::read_excel(filename, sheet=sheet_num, 
                                        col_names = FALSE) )
    ROWS = plate_dimensions_from_wells(number_of_wells)$rows
    COLS = plate_dimensions_from_wells(number_of_wells)$cols
    idx_conc = which(  tolower( xl[,1] ) == concentration_identifier )
    idx_treat = which( tolower( xl[,1] ) == treatment_identifier )
    idx_line = which(  tolower( xl[,1] ) == sample_identifier )
    if(length(idx_conc)==0){
       stop(paste("parameter concentration_identifier",concentration_identifier,
                   "not present in first column of Excel file",filename))
    }
    if(length(idx_treat)==0){
        stop(paste("parameter treatment_identifier",treatment_identifier,
                   "not present in first column of Excel file",filename))
    }
    if(length(idx_line)==0){
        stop(paste("parameter sample_identifier",sample_identifier,
                   "not present in first column of Excel file",filename))
    }
    if( dim(xl)[2] < COLS+1 ){
        stop( paste("fewer columns in plate map than expected, there should be",
                    COLS, "columns for a", number_of_wells,"well plate") )
    }
    
    rownames_conc =  trimws(xl[ (idx_conc+1) : (idx_conc+ROWS),1])
    rownames_treat = trimws(xl[ (idx_treat+1) : (idx_treat+ROWS),1])
    rownames_line = trimws(xl[ (idx_line+1) : (idx_line+ROWS),1])
    # asking for a value out of bounds for what read_excel() thinks of as the 
    # spreadsheet's bounds will return a NA
    if( sum(is.na(rownames_conc))>0 | 
        length(rownames_conc) != ROWS | 
        sum(rownames_conc != abc[1:ROWS] ) != 0){
        stop(paste("Expecting but did not see",ROWS,"row identifiers with",
                  "values A through", abc[ROWS], "for concentration plate map"))
    }
    
    if( sum( is.na(rownames_treat))>0 | 
        length(rownames_treat) != ROWS | 
        sum(rownames_treat != abc[1:ROWS] ) != 0 ){
        stop(paste("Expecting but did not see",ROWS,"row identifiers with",
                   "values A through", abc[ROWS], "for treatment plate map"))
    }    
    
    if( sum( is.na(rownames_line)>0) | 
        length(rownames_line) != ROWS | 
        sum( rownames_line != abc[1:ROWS] ) != 0 ){
        stop(paste("Expecting but did not see",ROWS,"row identifiers with",
                   "values A through", abc[ROWS], "for cell line plate map"))
    }    
    map_conc = data.frame(  xl[(idx_conc+1) : (idx_conc+ROWS),   2:(COLS+1)] )
    map_treat = data.frame( xl[(idx_treat+1) : (idx_treat+ROWS), 2:(COLS+1)] )
    map_line = data.frame(  xl[(idx_line+1) : (idx_line+ROWS),   2:(COLS+1)] )
    for( i in length(na.strings)){
        map_conc[ map_conc==na.strings[i] ] = NA   
        map_line[ map_line==na.strings[i] ] = NA
        map_treat[ map_treat==na.strings[i] ] = NA
    }
    map_conc = data.frame( data.matrix(map_conc) )
    names(map_conc) = 1:COLS
    names(map_treat) = 1:COLS
    names(map_line) = 1:COLS
    rownames(map_conc) = abc[1:ROWS]
    rownames(map_treat) = abc[1:ROWS]
    rownames(map_line) = abc[1:ROWS]

    if( ! is.numeric(map_conc[!is.na(map_conc)]) )
        stop( paste("All fields for the concentration matrix must be either a",
                    "number or the value NA")) 
    list(concentration=map_conc, treatment=map_treat, sample_type=map_line)
}

write.output = function(fn_out, o){
    # write values in vector o to fn_out, one per line
    for(i in 1:length(o) ){
        if(i==1)
            write(o[i], file=fn_out)
        else
            write(o[i], file=fn_out, append=TRUE)
    }
}

#' Write a loaded plate map object in XML format
#' 
#' @param fn A complete path to the file to write.
#' @param plate_map A list with matrixes labeled concentration, treatment, and 
#' sample_type.
#' @param passage Number for cell passage.
#' @param seeding_density Number for seeding density.
#' @param density_units Text string indicating the density units.
#' @param concentration_units Text string indicating the treatment concentration 
#' units.
#' @return none
#' @examples
#' pkg = "HTDoseResponseCurve"
#' fn_map = system.file("extdata", "sample_data_96_platemap.xlsx",package=pkg)
#' plate_map = read_platemap_from_excel( fn_map, number_of_wells=96 )
#' #write_platemap_to_XML( fn_map, plate_map, passage=1, seeding_density=1000,
#' #  density_units="thousands", concentration_units="ng/ul")
#' @export
write_platemap_to_XML = function( fn, plate_map, passage, seeding_density, 
                                  density_units, concentration_units){
    o = c()
    o=c(o,'<?xml version="1.0" encoding="utf-8"?>')
    o=c(o,'<plateMap>')
    o=c(o,'  <referenceItemManager>')
    o=c(o,'    <referenceItems>' )
    treatments = get_treatments(plate_map)
    for(i in 1:length(treatments)){
        tt = treatments[i]
        if( tt != "" ){
            x=paste('      <referenceItem type="Compound" description="',tt,
                    '" displayName="', tt,'" colorArgb="-12523200" />',
                    collapse='',sep='')
            o=c(o,x)
        }
    }
    sample_types = get_sample_types(plate_map)
    for(i in 1:length(sample_types)){
        st = sample_types[i]
        if( st!="" ){
            x=paste('      <referenceItem type="CellType" description="',st,
                    '" displayName="', st,'" colorArgb="-12523200" />',
                    collapse='',sep='')
            o=c(o,x)
        }
    }
    o=c(o,'    </referenceItems>')
    o=c(o,'  </referenceItemManager>')    
    
    N_ROW = dim(plate_map$treatment)[1]
    N_COL = dim(plate_map$treatment)[2]
    
    o=c(o,'<wellStore>')
    o=c(o,'  <wells>')
    for(rr in 1:N_ROW){
        for(cc in 1:N_COL){
            cur_t = plate_map$treatment[rr,cc]
            cur_s = plate_map$sample_type[rr,cc]
            cur_c = plate_map$concentration[rr,cc]
            if( cur_t==""){
                x = paste('      <well row="',rr-1,'" col="', cc-1, '" />',
                          collapse="", sep="") 
                o = c(o,x)
            }
            else{
                o = c(o, paste('      <well row="',rr-1,'" col="',cc-1, '" >',
                               collapse="", sep="") )
              o=c(o,      '        <items>')
              o=c(o,paste('          <wellItem type="Compound" concentration="',
                 cur_c,'" concentrationUnits="',concentration_units,'">',sep='', 
                            collapse=''))
              o=c(o,paste('            <referenceItem type="Compound" ',
                          'description="', cur_t,'" displayName="',cur_t,
                          '" colorArgb="-12523200" />', collapse='', sep=''))
                o = c(o,      '          </wellItem>')
                o = c(o,paste('          <wellItem type="CellType" passage="', 
                              passage, '" seedingDensity="', seeding_density, 
                              '" seedingDensityUnits="', density_units, '">', 
                              sep='', collapse=''))
                o = c(o,paste('            <referenceItem type="CellType" ', 
                              'description="', cur_s,'" displayName="',cur_s,
                              '" colorArgb="-1490880" />', 
                              collapse='', sep=''))
                o = c(o,      '          </wellItem>')
                o = c(o,      '        </items>')
                o = c(o,      '      </well>')
            }       
        }
    }
    o=c(o,'    </wells>')
    o=c(o,'  </wellStore>')
    o=c(o,'</plateMap>')
    
    write.output(fn, o)
}


#' Read Incucyte-formatted Excel document for a plate at one or more time points
#'
#' \code{read_plates_from_Incucyte_export} loads data that have been exported to 
#' Excel from the Incucyte platform and returns a data frame of the results, 
#' with additional columns for the timestamp, hour, and a user-defined plate 
#' identifier. Each call to this function will load the exported values from one 
#' or more time points for a single plate map. Individual plate maps can be 
#' distinguished in later analysis by the user-provided plate_id parameter.
#' 
#' The number of columns will be the columns in the plate plus three 
#' additional columns: timestamp, plate_id, and hours. The number of rows will 
#' be the number of rows in the plate times the number of plates in the Excel 
#' sheet. Excel sheets are read using Hadley Wickham's readxl package.
#' 
#' If the Incucyte instrument only exports a partial plate (e.g. if it does not 
#' export every row) this function will including padding for rows and columns 
#' that are not in the extract.
#' 
#' @param path_to_file Path to Excel file generated by the Incucyte software
#' @param plate_id A user-defined string to identify this plate. A plate can be 
#' imaged at more than one time point.
#' @param number_of_wells The number of wells per plate, must be 384 or 96
#' @param sheet_num The index of the Excel sheet to read. Defaults to 1.
#' @param plate_mask Indicates which wells to read. Defaults to NA, read all 
#' wells. If a matrix of boolean values is passed, read only the wells where 
#' plate_mask is TRUE. The matrix must have the same dimensions as the wells 
#' being read.
#' @return A data frame where columns are data, plate_id, hour. 
#' @examples 
#' pkg = "HTDoseResponseCurve"
#' fn_data_Excel = system.file("extdata", "sample_data_384.xlsx", package = pkg)
#' plate_data = read_plates_from_Incucyte_export( fn_data_Excel, "p1", 
#'                                                number_of_wells=384)
#' @export
read_plates_from_Incucyte_export = function( path_to_file, plate_id, 
                                             number_of_wells, sheet_num=1, 
                                             plate_mask = NA){
    
    expected_rows = plate_dimensions_from_wells(number_of_wells)$rows
    expected_cols = plate_dimensions_from_wells(number_of_wells)$cols
    if( is.na(plate_mask) ){
        plate_mask = matrix(TRUE, nrow= expected_rows, ncol=expected_cols )
    }else{
        if( dim(plate_mask)[1] != expected_rows |
            dim(plate_mask)[2] != expected_cols ){
            stop("plate_mask dimension different from data plate dimension")   
        }
    }
    
    xl = data.frame( readxl::read_excel(path_to_file, sheet=sheet_num) )
    n_xl_cols = dim(xl)[2]
    # incucytye can't be trusted to write all rows or all columns
    rows_timestamp = which( xl[,1] == "Time Stamp:" )
    rows_a1 = which( xl[,1]=="A")
    if( length(rows_timestamp)>1 ){
        n_plate_data_rows = rows_timestamp[2]-rows_a1[1]-1
    }else{
        n_plate_data_rows = dim(xl)[1]-rows_a1[1]+1
    }
    highest_col_letter = xl[ rows_a1[1]:dim(xl)[1],1]
    for(i in 1:length(rows_a1)){
        hour = as.numeric( xl[rows_timestamp[i], 4] )
        plate = create_empty_plate( number_of_wells, hour, plate_id)
        col_a1 = 2
        row_last = rows_a1[i] + n_plate_data_rows
        col_last = n_xl_cols
        plate_read = data.matrix(xl[rows_a1[i] : row_last, col_a1 : col_last])
        plate[1:dim(plate_read)[1], 1:dim(plate_read)[2]] = plate_read
        vals = plate[1:expected_rows, 1:expected_cols]
        vals[!plate_mask]  = NA
        plate[1:expected_rows, 1:expected_cols] = vals
        if(i==1){
            plates=plate
        }
        else{
            plates = rbind(plates, plate)
        }
    }
    plates
}


#' Read a dataset from a text file
#' 
#' The text file may be separated by commas, tabs, or any other commonly used 
#' delimiter. The first row must have the following elements: 
#' \itemize{
#'  \item{sample_type}
#'  \item{treatment}
#'  \item{concentration}
#'  \item{value}
#' }
#' 
#' The text file may optionally contain a column called "hour", indicating the 
#' time point of the measurement. If "hour" is present the value must be a 
#' number.
#' 
#' @param filepath Path to the file to load. If the full path is not given, R 
#' will assume the file is in the current working directory.
#' @param negative_control A designation for the negative controls in this 
#' dataset, if they exist. Value may be NA, a number, a string, or a data frame.
#' 
#' \itemize{
#'  \item{NA: Use when there are no negative control measurements. The contents 
#'    of the column named 'value_normalized' will be copied from the contents of 
#'    the column named 'value'. }
#'  \item{Number: Use when each treatment has been labeled with a concentration 
#'    (typically 0) that indicates the vehicle control. Each treatment must 
#'    contain one or more observations with this concentration, and these 
#'    observations will be the negative controls.}
#'  \item{string: Use when a single set of observations is a universal control. 
#'    The treatment whose name matches the string is the universal 
#'    negative control all of the data.}
#'  \item{data frame: Use when more than one negative control exists, and you 
#'    have to map different treatments to a particular negative control. The 
#'    data frame must have names 'drug' and 'vehicle', and the data frame will 
#'    map match treatments in the 'drug' column to those in the 
#'    'vehicle' column.}
#' }
#' @param plate_id A string identifying this set of observations. Defaults to 
#' "plate_1".
#' @param sep A character that separates elements in the text file. Fefaults to 
#' tab.
#' @export
#' @return A data frame for the dataset matching the output of the 
#' \code{\link{create_dataset}} function.
#' @examples 
#' pkg = "HTDoseResponseCurve"
#' fn_txt = system.file("extdata", "sample_data_1.txt", package = pkg)
#' ds = read_dataset( fn_txt, negative_control="DMSO" )
#' @importFrom utils read.table
#' @seealso \code{\link{create_dataset}} 
read_dataset = function(filepath, negative_control=NA, plate_id="plate_1", 
                        sep="\t"){
    
    ds = utils::read.table(filepath, sep=sep,header=TRUE,stringsAsFactors=FALSE)
    if( ! "sample_type" %in% names(ds) ){
        stop( "dataset file must contain a column labeled 'sample_type'")
    }
    if( ! "treatment" %in% names(ds) ){
        stop( "dataset file must contain a column labeled 'treatment'")
    }
    if( ! "concentration" %in% names(ds) ){
        stop( "dataset file must contain a column labeled 'concentration'")
    }
    if( ! "value" %in% names(ds) ){
        stop( "dataset file must contain a column labeled 'value'")
    }
    hours=0
    if( ! "hour" %in% names(ds) ){
        ds$hour = rep(0, dim(ds)[1])   
    }
    hours = ds$hour
    if( sum(! is.numeric(ds$hour) )>0 ){
        stop("'hour' column must be numeric")
    }
    if( "concentration_2" %in% names(ds) ){
        concentrations_2 = ds$concentration_2   
    }else{
        concentrations_2 = NULL
    }
    if( "treatment_2" %in% names(ds) ){
        treatments_2 = ds$treatment_2   
    }else{
        treatments_2 = NULL
    }    
    create_dataset(sample_types=ds$sample_type, 
                   treatments = ds$treatment,
                   treatments_2 = treatments_2,
                   concentrations = ds$concentration, 
                   concentrations_2 = concentrations_2,
                   hours = hours,
                   values = ds$value,
                   plate_id= plate_id,
                   negative_control = negative_control)
}