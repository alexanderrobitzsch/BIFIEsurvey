## File Name: load_BIFIEdata_files_cat_print_file_name.R
## File Version: 0.02


load_BIFIEdata_files_cat_print_file_name <- function( file )
{
    cat(paste0( "- Read ", file, "\n") )
    utils::flush.console()
}
