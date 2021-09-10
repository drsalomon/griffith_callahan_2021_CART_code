"""
=================================================================================================
Kenneth P. Callahan

17 July 2021

=================================================================================================
Python >= 3.8.5

extract_kegg_go_data.py

This scripts can take the same input as the pepdep_to_gct.py, but rather than providing input
for PTM-SEA, this script parses the elements of specific columns of the input file and
makes bar plots using them. Currently, the columns to look for are:

   Renamed Col       Meaning
---------------------------------------------------------------------
kegg            <-> Kyoto Encyclopedia of Genes and Genomes
GO_molecular    <-> Gene Ontology Molecular Function Annotaitons
GO_localization <-> Gene Ontology Localization Annotations
---------------------------------------------------------------------

Note that the renamed columns are after renaming the PeptideDepot Data Dump columns.

=================================================================================================
Arguments:

args[1] (data_file)        -> The path to a PeptideDepot data dump (excel file).
args[2] (metadata_file)    -> The path to the sample metadata file, which contains the
                              required columns (reordered), renamed columns (same order as
                              required columns), peakarea_manual_1/2 conditions, timepoints/conditions
args[3] (output_dir)       -> The directory to write all outputs to. This directory is
                              not required to exist prior to running the script.
args[4] (sample_delim)     -> The delimiter for the sample metadata text file. Currently, I am using
                              a comma (,).
args[5] (ascore_value)     -> The minimum ambiguity score (A-score) to be considered in PTM-SEA

=================================================================================================
Outputs:


output_directory
    |
    |
    |------log
    |       |
    |       |------<condition_1>
    |       |           |
    |       |           |-----------<pvalue>  -> GO/KEGG Barcharts
    |       |           |-----------<qvalue>  -> GO/KEGG Barcharts
    |       |
    |       |------<condition_2>
    |                   |
    |                   |-----------<pvalue>  -> GO/KEGG Barcharts
    |                   |-----------<qvalue>  -> GO/KEGG Barcharts
    |
    |------raw
            |
            |------<condition_1>
            |           |
            |           |-----------<pvalue>  -> GO/KEGG Barcharts
            |           |-----------<qvalue>  -> GO/KEGG Barcharts
            |
            |------<condition_2>
                        |
                        |-----------<pvalue>  -> GO/KEGG Barcharts
                        |-----------<qvalue>  -> GO/KEGG Barcharts

*** Note that the conditions will depend on which columns you choose to keep, ***
*** and how you chose to define the conditions in the sample metadata file.   ***

=================================================================================================
"""
############################################################################################################
#
#     Importables

# peptide_to_gct is a script used for formatting PeptideDepot data dumps
# in such a way that is compatible with PTM-SEA. This script borrows a lot
# from that script, as it uses the same input file and parameters.
from py_scripts import pepdep_to_gct as ptg

# General helpers has a number of functions I use frequently in
# my scripts. They are all placed in that module purely for
# convenience and generalizability.
import py_scripts.helpers.general_helpers as gh

# The Pandas Helper file has scripts that help manage Pandas
# DataFrames, and perform various actions on lists of dataframes
import py_scripts.helpers.pandas_helpers as ph

# MPL Plotting Helpers has plotting functions, in particular the
# bar plots.
import py_scripts.helpers.mpl_plotting_helpers as mph

# The random module is only used if using random selections
import random

# os is used in order to perform operating system level tasks.
import os
import sys


#
#
############################################################################################################
#
#     Global Variables

# The field_strings list determines which values to keep
# from the ttest output. Currently, we are keeping both
# the pvalue column and the qvalue column
field_strings = ["pvalue","qvalue"]

# This dictionary holds the settings to be applied for
# T-Tests. For more information on what this dictionary
# controls, see homebrew_stats.pairwise_t() and
# homebrew_stats.ttest()
df_t_settings = {"identity_column"  : "peptide",
                 "comparisons"      : "all",
                 "pw_ttest_optargs" : {"t_type" : "student"},
                 "qvalue"           : True,
                 "storey_pi0"       : 1}

# This list holds the strings representing the columns of
# interest for creating bar charts.
barplot_col_strings = ["Kegg unique index",
                       "Go mol process unique index",
                       "Go location unique index",
                       "Go bio process unqiue index"]

# This is a list of all column headers from PeptideDepot
# that can be ignored for this program. If they are
# found in the sample metadata, they will simply be ignored.
ignore_cols = ["all protein name index",
               "peptide sequence GCT format centered on 1st site",
               "peptide sequence GCT format centered on 2nd site",
               "peptide sequence GCT format centered on 3rd site",
#               "Kegg unique index",
#               "Go mol process unique index",
               "Go bio process unique index",
#               "Go location unique index",
               "swissprot accession",
               "string accession",
               "hprd accession",
               "UNIPROT Gene Name",
               "UNIPROT accession number",
               "gi species adjust",
               "xcorr max across timepoints",
               "qvalues for SILAC timepoint1",
               "qvalues for SILAC timepoint2",
               "qvalues for SILAC timepoint3",
               "heatmap qvalues for labelfree 1 timepoint1",
               "heatmap qvalues for labelfree 1 timepoint2",
               "heatmap qvalues for labelfree 1 timepoint3",
               "heatmap qvalues for labelfree 2 timepoint1",
               "heatmap qvalues for labelfree 2 timepoint2",
               "heatmap qvalues for labelfree 2 timepoint3",
#               "protein name manual",
               "assigned sequence phospho probabilities highlight all sites",
               "charge state peptide",
               "phosphosite annotated",
               "accession number for psite"]

#
#
############################################################################################################
#
#     Functions

def format_strcols(*lists,
                   delim = "\n"):
    """
    =================================================================================================
    format_strcols(*lists, delim)
    
    This function is meant to take an arbitrary number of lists that came from the column of an
    excel file and return a list of all strings from the column.
    
    =================================================================================================
    Arguments:
    
    lists  ->  An arbitrary number of lists of strings.
    delim  ->  A string describing the character separating elements in the column.
    
    =================================================================================================
    Returns: A list of strings
    
    =================================================================================================
    """
    # Initialize the columns list
    cols = []
    # Loop over each list in the input lists
    for a_list in lists:
        # First, filter out all nan values in the list
        a_list = [item for item in a_list if item == item and item.lower() != "nan"]
        # Then, split each element on the delimiter
        a_list = [item.split(delim) for item in a_list]
        # and unpack the list, and add that list to the cols
        cols.append(gh.unpack_list(a_list))
    # At the end, return the cols list.
    return cols

def format_counts(*dicts):
    """
    =================================================================================================
    format_counts(*dicts)
    
    This function is meant to take an arbitrary number of dictionaries created using gh.count()
    and merge them together, with the keys from the dicts and the values as lists of counts.
    
    =================================================================================================
    Arguments:
    
    dicts  ->  An arbitrary number of counting dictionaries
    
    =================================================================================================
    Returns: The "xlabels" (a list of the keys from the dictionaries) and the "yvals" (a list of
             lists of values that come from the counts in the input dictionary. Sublists are
             index paired with the input dictionaries).
    
    =================================================================================================
    """
    # Merge the count dictionaries on the keys
    counts = gh.merge_dicts(*dicts)
    # Turn the dictionary into a list of tuples, with the key value pairs.
    count_tups = [(key, value) for key, value in counts.items()]
    # Then, sort this list based on the maximum value in each count list.
    # This is to ensure the bars are in (mostly) descending order.
    sort_count = sorted(count_tups, key = lambda x: x[1][x[1].index(max(x[1]))], reverse= True)
    # The labels for each group are in the zeroeth element of each tuple
    xlabels = [category[0] for category in sort_count]
    # The y-values are in the first element of each tuple. However, matplotlib
    # takes the transpose of specific points, so we transpose this list. The
    # resulting list will contain sublist where sublist i corresponds to
    # dictionary i, where missing values have been replaced by zeroes.
    yvals = gh.transpose(*[category[1] for category in sort_count])
    # At the end, return both the xlabels and the yvalues.
    return xlabels, yvals

def merge_fieldvals_and_main(main_dict,
                             stats_dict,
                             main_df):
    """
    =================================================================================================
    merge_fieldvals_and_main(main_dict, stats_dict, main_df)
    
    This function takes a dictionary of dataframes, a dictionary of statistical outputs, and the
    original dataframe used for the construction of both dictionaries and returns a dictionary
    where the statistical tests in stats_dict are merged with the main_df.
    
    This function assumes that the main dataframe and the dataframes from stats dict have the
    same indices. It is suggested that the index of main_df is set, inplace, to the identity column.
    
    =================================================================================================
    Arguments:
    
    main_dict    ->  A dictionary containing dataframes where all values besides the ID column and
                     data columns have been removed. The keys in this dataframe should also be a
                     key in the stats_dict.
    stats_dict   ->  A dictionary containing dictionaries containing statistical tests on the
                     dataframes from main_dict.
    main_df      ->  A pandas DataFrame that is the unmodified dataframe used for both main_dict
                     and stats_dict.
    
    =================================================================================================
    Returns: A new dictionary where the statistical tests from stats_dict are merged index-wise
             with the main dict.
    
    =================================================================================================
    """
    # Initialize the newdict variable, which will hold the outputs.
    newdict = {}
    # Loop over the keys and values in the main dictionary.
    for key_main, val_main in main_dict.items():
        # Initialize a subdictioary in newdict under the key_main
        newdict[key_main] = {}
        # And loop over the keys and values of the subdictionary in stats_dict
        # under this main key.
        for key_sub, val_sub in stats_dict[key_main].items():
            # Use ph.df_combine() to combine the main_df with the value df
            newdict[key_main][key_sub] = ph.df_combine(main_df, val_sub)
    # And once all combinations are complete, return the newdict.
    return newdict

def extract_heads(a_stats_dict):
    """
    =================================================================================================
    extract_heads(a_stats_dict)
    
    This function is meant to extract the column headers from dataframes in a statistics dictionary.
    The dataframes in the statistics dictionary are asumed to have undergone ph.df_parse_fieldvalues,
    which takes the comparison groups and the statistics associated with them and makes columns
    based on the groups, the statistics values, and the index of the dataframe.
    
    =================================================================================================
    Arguments:
    
    a_stats_dict  ->  A dictionary of dataframes which represent the results of statistical tests
                      on the indices.
    
    =================================================================================================
    Returns: A list of strings representing the column headers from the dictionary.
    
    =================================================================================================
    """
    # Loop over the keys and values in the statistics dictionary
    for key, value in a_stats_dict.items():
        # if the value is a dictionary
        if type(value) == dict:
            # then run this fucntion recursively
            return extract_heads(value)
        # Otherwise
        else:
            # return a list of the column headers.
            return list(value.columns.values)
        
def count_proteins(a_dict_of_dfs,
                   comp_cols,
                   sample_metadata_file,
                   meta_delim = ",",
                   outpath = "",
                   threshold = 0.05):
    """
    =================================================================================================
    count_proteins(*args, **kwargs)
    
    The setup here is almost wrap_barplots() verbatum
    
    =================================================================================================
    Arguments:
    
    a_dict_of_dfs         ->  A dictionary of dataframes, which statistics have been added to.
    a_list_of_cols        ->  A list of the column names that we would like to turn into barplots.
    comp_cols             ->  A list of strings representing the comparisons represented by the bars
    outpath               ->  A string represnting the output filepath for the image
    threshold             ->  A float between 0 and 1 representing the maximum value for counting
                              consideration in the barplot.
    sample_metadata_file  ->  The name of the sample metadata file
    meta_delim            ->  The delimiter used in that file, default is comma
    
    =================================================================================================
    Returns: None
    
    =================================================================================================
    """
    # Find the protein name manual column using the sample metadata
    prot_col = ptg.find_col(sample_metadata_file,
                            "protein name manual",   # Should get a list of all possible values for this
                            delim = meta_delim)
    # If no protein name manual was found
    if prot_col.lower() == "none":
        # Then return None, as there is no protein names to count
        return None
    # Loop over the keys and values in the input dictionary
    for key, value in a_dict_of_dfs.items():
        if type(value) == dict:
            # If the value is a dictionary, then run this function
            # recursively using the value, and update the path
            count_proteins(value,
                           comp_cols,
                           sample_metadata_file,
                           meta_delim = meta_delim,
                           outpath = os.path.join(outpath, key),
                           threshold = threshold)
        # If the current value is not a dictionary, then we
        # assume it is a dataframe.
        else:
            # If the output path is not currently in existance
            if not os.path.exists(os.path.join(outpath, key)):
                # Then make the directory so we can write files to it
                os.makedirs(os.path.join(outpath,key))
            # Otherwise,loop over the comparison labels. We want to see
            # how many significantly changing proteins exist per comparison
            for comp in comp_cols:
                # First, filter the current dataframe based on a threshold value
                filtered_val = value[value[comp].astype(float) < threshold]
                # Then count the names in the protein column
                counted = gh.count(list(filtered_val[prot_col]))
                # And turn that dictionary into a list of lists
                counted = [[key,value] for key, value in counted.items()]
                # Next, add the column headers, and sort the list based on the names
                counted = [["protein name", "count"]] + sorted(counted, key = lambda x: x[0])
                # and convert all sublists into strings, tab separated
                counted = [gh.list_to_str(item,
                                          delimiter = "\t",
                                          newline = True) for item in counted]
                # Then write the output file with the counted proteins.
                gh.write_outfile(counted,
                                 os.path.join(outpath,key,f"protein_counts_{comp}_sigthresh.txt"))
    return None

def wrap_barplots(a_dict_of_dfs,
                  a_list_of_cols,
                  comp_cols,
                  outpath = "",
                  threshold = 0.05,
                  max_bars = 25,
                  colour_list = ["blues", "pinks", "greens", "blacks", "oranges", "all"],
                  index = 0,
                  increment = True):
    """
    =================================================================================================
    wrap_barplots(*args, **kwargs)
    
    This function is meant to create a number of barplots at the same time.
    
    =================================================================================================
    Arguments:
    
    a_dict_of_dfs   ->  A dictionary of dataframes, which statistics have been added to.
    a_list_of_cols  ->  A list of the column names that we would like to turn into barplots.
    comp_cols       ->  A list of strings representing the comparisons represented by the bars
    outpath         ->  A string represnting the output filepath for the image
    threshold       ->  A float between 0 and 1 representing the maximum value for counting
                        consideration in the barplot.
    max_bars        ->  An integer representing the maximum number of categories permitted on
                        the plot.
    colour_list     ->  A list of colour group identities. For more info, see the colours
                        dictionary in mpl_plotting_helpers
    index           ->  The value for the current index of the loop.
    increment       ->  A boolean determining whether or not to increment the index.
    
    =================================================================================================
    Returns: None
    
    =================================================================================================
    """
    # Initialize the barplot_settings dictionary. For more information on
    # the acceptable settings for barplots, see mpl_plotting_helpers.plot_bars()
    barplot_settings = {"col_labels" : comp_cols,
                        "colour_choice" : "centered",
                        "separation" : 0.3,
                        "img_type" : "pdf",
                        "show" : False,
                        "subplot_args" : {"figsize" : (24,12)}}
    # Loop over the keys and values in the input dictionary
    for key, value in a_dict_of_dfs.items():
        # If the type of value is a dictionary,
        if type(value) == dict:
            # Then try to run this function recursively using the value.
            # Also, change the increment to False, we want the colours
            # within groups to stay the same.
            wrap_barplots(value,
                          a_list_of_cols,
                          comp_cols,
                          outpath = os.path.join(outpath, key),
                          threshold = threshold,
                          max_bars = max_bars,
                          colour_list = colour_list,
                          index = index,
                          increment = False)
        # If the value is not a dictionary, then we assume
        # it is a dataframe.
        else:
            # Set the random seed to index plus 900. If we are choosing
            # random colours from 'all', we want the colours to be consistent
            # within a group.
            random.seed(900 + index)
            # Next, set the colour_type of barplot_settings to the
            # colour_list at index
            barplot_settings["colour_type"] = colour_list[index]
            # Next, loop over the list of desired columns. These are the columns
            # we would like to make a barplot for.
            for col in a_list_of_cols:
                # Check to see if the output filepath exists or not.
                if not os.path.exists(os.path.join(outpath,key)):
                    # If not, create the path
                    os.makedirs(os.path.join(outpath,key))
                # Now that the output path exists, we can update the barplot settings
                # Make the image name the column name with the barplot tag
                barplot_settings["img_name"] = os.path.join(outpath, key, f"{col}_barplot")
                # And update hte set_kwargs with the proper values for x,y and title.
                barplot_settings["set_kwargs"] = {"xlabel" : f"{col.upper()} Categories",
                                                  "ylabel" : f"Counts",
                                                  "title"  : f"{col.upper()} Representation for {key} < {threshold}"}
                # Next, filter the current dataframe based on the threshold value
                # for each column in the comparisons strings list.
                all_filts = [value[value[comp].astype(float) < threshold] for comp in comp_cols]
                # Next, turn the column containing the barplot values for
                # plotting into lists, for each filtered dataframe
                all_vals = [list(filt_df[col]) for filt_df in all_filts]
                # Then, format the strings in each of these lists
                all_vals = format_strcols(*all_vals)
                # and count all of the lists.
                ind_counts = gh.count_lists(*all_vals)
                # Then, format the count dictionaries into xlabels and yvalues
                all_counts = format_counts(*ind_counts)
                # And plot the barplot for this particular column.
                mph.plot_bars(all_counts[0][:max_bars],
                              [[item[i] for i in range(len(item)) if i < max_bars]
                              for item in all_counts[1]],
                              **barplot_settings)
        # If the increment boolean is True
        if increment:
            # Then increment the index.
            index += 1
    # Return None once the function is completed.
    return None

#
#
############################################################################################################
#
#     main() functions

def imp_main(pepdep_xls_file,
             sample_metadata_file,
             meta_delim = ",",
             ascore = 13):
    """
    =================================================================================================
    imp_main(pepdep_xls_file, sample_metadata_file, meta_delim, ascore)
    
    This function is the IMPortable MAIN function, and as it's name suggests, can be used as when
    imported. This function will provide barplots for three particular columns of a PeptideDepot
    data dump:
    
    KEGG
    GO molecular processes
    GO Localization
    
    The names may be a little bit wonky, so the output is an editable PDF file.
    
    =================================================================================================
    Arguments:
    
    args[1] (data_file)        -> The path to a PeptideDepot data dump (excel file).
    args[2] (metadata_file)    -> The path to the sample metadata file, which contains the
                                  required columns (reordered), renamed columns (same order as
                                  required columns), peakarea_manual_1/2 conditions, timepoints/conditions
    args[3] (sample_delim)     -> The delimiter for the sample metadata text file. Currently, I am using
                                  a comma (,).
    args[4] (ascore_value)     -> The minimum ambiguity score (A-score) to be considered in PTM-SEA
    
    =================================================================================================
    Returns: None
    
    =================================================================================================
    """
    # State which variables should be used from the global scope.
    global ignore_cols         # Columns to ignore from peptide depot xls file
    global field_strings       # List of strings for 'field'-values
    global barplot_col_strings # List of headers for barplot columns
    global df_t_settings       # Settings for df_ttests()
    # Use the path to the peptide_depot xls file as the output path
    path = os.path.abspath(os.path.dirname(pepdep_xls_file))
    # Check the arguments for validity
    print("Checking the inputs for validity....\n")
    print(pepdep_xls_file,
                             sample_metadata_file,
                             path,
                             meta_delim,
                             ascore)
    args = ptg.check_sysargs(pepdep_xls_file,
                             sample_metadata_file,
                             path,
                             meta_delim,
                             ascore,
                             ported = True)
    # Make a list of the columns of interest by finding which columns
    # from the metadata file correspond to the barplot columns
    print("Finding the non-data columns...\n")
    cols_of_interest = [ptg.find_col(args[1],
                                     col,
                                     delim = args[3]) for col in barplot_col_strings]
    # Do the same to find the ascore column
    ascore_col = [ptg.find_col(args[1],
                               "ascore max",            # Should get a list of all possible values for this
                               delim = args[3])]
    # and the column with the assigned sequence
    assigned_col = [ptg.find_col(args[1],
                                 "assigned sequence",   # Should get a list of all possible values for this
                                 delim = args[3])]
    prot_col = [ptg.find_col(args[1],"protein name manual", delim = args[3])]
    print("Columns found:\n")
    for item in cols_of_interest + ascore_col + assigned_col + prot_col:
        print(item)
    print("\n")
    # If the assigned column was not found
    if assigned_col == "none":
        # Raise a value error, as this column is required
        raise ValueError("The peptide sequence assigned cannot be found...")
    # Once the filepath is determined, make the new directory string for outputs
    path = os.path.join(args[2], "kegg_go_outfiles")
    print("Reading the input data file...")
    # Then read the peptide depot data dump, ignoring the columns listed globally
    pepdep_xls_df, experiment_metadata = ptg.read_pepdep_dump(args[0],
                                                              args[1],
                                                              meta_delim = args[3],
                                                              ignore = ignore_cols)
    # And filter the peptide depot data dump based on the given ascore
    print(f"Filtering the input data file by an ambiguity score of {args[4]}...")
    pepdep_xls_df = ptg.filt_xls_by_ascore(pepdep_xls_df,
                                           args[1],
                                           delim = args[3],
                                           ascore_val = args[4])
    # Then, find all columns that are not data columns
    id_cols = [item for item in list(pepdep_xls_df.columns.values) if not any([True for _ in experiment_metadata[0] if _ in item])]
    # and clean up the excel file using clean_pepdep_data, and return dataframes
    data_dict, parsing_strs = ptg.clean_pepdep_data(pepdep_xls_df,
                                                    experiment_metadata,
                                                    id_col = id_cols,
                                                    value_df = True)
    # Use dictionary comprehension to remove all non-id and non-data columns.
    data_dict_vals = {key : value[[col for col in list(value.columns.values)  if col not in ascore_col + cols_of_interest + prot_col]]
                     for key, value in list(data_dict.items())}
    # And perform statistics on those dataframes
    print("Performing statistics on the raw data...")
    stats_dict_vals = ph.df_ttest_on_dict(data_dict_vals,
                                          experiment_metadata[1],
                                          parsing_strs,
                                          fields = field_strings,
                                          ttest_settings = df_t_settings,
                                          filename = "statisitcs",
                                          path = os.path.join(path, "raw") )
    # Then, log transform the data
    data_dict_log = ptg.log_data(data_dict_vals,
                                 id_col = assigned_col)
    print("Performing statistics on the log data...")
    # and perform statistics on those log transformed dataframes
    stats_dict_log = ph.df_ttest_on_dict(data_dict_log,
                                         experiment_metadata[1],
                                         parsing_strs,
                                         fields = field_strings,
                                         ttest_settings = df_t_settings,
                                         filename = "statistics_log",
                                         path = os.path.join(path, "log"))
    # Then extract the headers from the stats dictionary dataframes,
    # these are the comparisons used for statistics
    comp_cols = extract_heads(stats_dict_log)
    # and reset the index of the main dataframe to the assigned column
    pepdep_xls_df.set_index(assigned_col[0], inplace = True)
    # Then, merge the raw data dataframes with the raw data statistics
    data_dict_vals = merge_fieldvals_and_main(data_dict_vals,
                                              stats_dict_vals,
                                              pepdep_xls_df)
    # And merge the log data dataframes with the log data statistics
    data_dict_log = merge_fieldvals_and_main(data_dict_log,
                                             stats_dict_log,
                                             pepdep_xls_df)
    print(f"Counting proteins based on threshold {0.05} for raw data...")
    count_proteins(data_dict_vals,
                   comp_cols,
                   args[1],
                   meta_delim = args[3],
                   outpath = os.path.join(path, "raw"))
    print(f"Counting proteins based on threshold {0.05} for log data...")
    count_proteins(data_dict_log,
                   comp_cols,
                   args[1],
                   meta_delim = args[3],
                   outpath = os.path.join(path, "log"))
    print(f"Plotting all bar plots to {args[2]}")
    # Then, use the dictionaries to make barplots.
    wrap_barplots(data_dict_vals,
                  cols_of_interest,
                  comp_cols,
                  outpath = os.path.join(path, "raw"))
    wrap_barplots(data_dict_log,
                  cols_of_interest,
                  comp_cols,
                  outpath = os.path.join(path, "log"))
    print("Complete!! :)\n")
    # and return None at the end.
    return None

def cmd_main():
    """
    =================================================================================================
    cmd_main()
    
    This function is the CoMmanD line MAIN function, and as it's name suggests, can be used when
    run from the command line. This function will provide barplots for three particular columns of
    a PeptideDepot data dump:
    
    KEGG
    GO molecular processes
    GO Localization
    
    The names may be a little bit wonky, so the output is an editable PDF file.
    
    =================================================================================================
    Arguments:
    
    args[1] (data_file)        -> The path to a PeptideDepot data dump (excel file).
    args[2] (metadata_file)    -> The path to the sample metadata file, which contains the
                                  required columns (reordered), renamed columns (same order as
                                  required columns), peakarea_manual_1/2 conditions, timepoints/conditions
    args[3] ()                 -> Any string. This will be immediately overwritten with the path to
                                  the peptide depot data dump.
    args[3] (sample_delim)     -> The delimiter for the sample metadata text file. Currently, I am using
                                  a comma (,).
    args[4] (ascore_value)     -> The minimum ambiguity score (A-score) to be considered in PTM-SEA
    
    =================================================================================================
    Returns: None
    
    =================================================================================================
    """
    # Get the command line arguments
    args = list(sys.argv)
    # and assign the third one to be the current working directory
    args[3] = os.getcwd()
    # and run check_sysargs from pepdep_to_gct.py
    print("Checking the inputs for validity....\n")
    args = ptg.check_sysargs(*args,
                             ported = False)
    # Then, reassign the third value to be the path to the peptide depot file.
    args[3] = os.path.abspath(os.path.dirname(args[1]))
    # Once all of the arguments have been checked,
    # determine which variables are defined globally
    global ignore_cols
    global field_strings
    global barplot_col_strings
    global df_t_settings
    # Make a list of the columns of interest by finding which columns
    # from the metadata file correspond to the barplot columns
    print("Finding the non-data columns...\n")
    cols_of_interest = [ptg.find_col(args[2],
                                     col,
                                     delim = args[4]) for col in barplot_col_strings]
    # Do the same to find the ascore column
    ascore_col = [ptg.find_col(args[2],
                               "ascore max",            # Should get a list of all possible values for this
                               delim = args[4])]
    # and the column with the assigned sequence
    assigned_col = [ptg.find_col(args[2],
                                 "assigned sequence",   # Should get a list of all possible values for this
                                 delim = args[4])]
    prot_col = [ptg.find_col(args[2],"protein name manual", delim = args[4])]
    print("Columns found:\n")
    for item in cols_of_interest + ascore_col + assigned_col + prot_col:
        print(item)
    print("\n")
    # If the assigned column was not found
    if assigned_col == "none":
        # Raise a value error, as this column is required
        raise ValueError("The peptide sequence assigned cannot be found...")
    # Once the filepath is determined, make the new directory string for outputs
    path = os.path.join(args[3], "kegg_go_outfiles")
    # Then read the peptide depot data dump, ignoring the columns listed globally
    print("Reading the input data file...")
    pepdep_xls_df, experiment_metadata = ptg.read_pepdep_dump(args[1],
                                                              args[2],
                                                              meta_delim = args[4],
                                                              ignore = ignore_cols)
    # And filter the peptide depot data dump based on the given ascore
    print(f"Filtering the input data file by an ambiguity score of {args[4]}...")
    pepdep_xls_df = ptg.filt_xls_by_ascore(pepdep_xls_df,
                                           args[2],
                                           delim = args[4],
                                           ascore_val = args[5])
    # Then, find all columns that are not data columns
    id_cols = [item for item in list(pepdep_xls_df.columns.values) if not any([True for _ in experiment_metadata[0] if _ in item])]
    # and clean up the excel file using clean_pepdep_data, and return dataframes
    data_dict, parsing_strs = ptg.clean_pepdep_data(pepdep_xls_df,
                                                    experiment_metadata,
                                                    id_col = id_cols,
                                                    value_df = True)
    # Use dictionary comprehension to remove all non-id and non-data columns.
    data_dict_vals = {key : value[[col for col in list(value.columns.values)  if col not in ascore_col + cols_of_interest+prot_col]]
                     for key, value in list(data_dict.items())}
    # And perform statistics on those dataframes
    print("Performing statistics on the raw data...")
    stats_dict_vals = ph.df_ttest_on_dict(data_dict_vals,
                                          experiment_metadata[1],
                                          parsing_strs,
                                          fields = field_strings,
                                          ttest_settings = df_t_settings,
                                          filename = "statisitcs",
                                          path = os.path.join(path, "raw") )
    # Then, log transform the data
    data_dict_log = ptg.log_data(data_dict_vals,
                                 id_col = assigned_col)
    # and perform statistics on those log transformed dataframes
    print("Performing statistics on the log data...")
    stats_dict_log = ph.df_ttest_on_dict(data_dict_log,
                                         experiment_metadata[1],
                                         parsing_strs,
                                         fields = field_strings,
                                         ttest_settings = df_t_settings,
                                         filename = "statistics_log",
                                         path = os.path.join(path, "log"))
    # Then extract the headers from the stats dictionary dataframes,
    # these are the comparisons used for statistics
    comp_cols = extract_heads(stats_dict_log)
    # and reset the index of the main dataframe to the assigned column
    pepdep_xls_df.set_index(assigned_col[0], inplace = True)
    # Then, merge the raw data dataframes with the raw data statistics
    data_dict_vals = merge_fieldvals_and_main(data_dict_vals,
                                              stats_dict_vals,
                                              pepdep_xls_df)
    # And merge the log data dataframes with the log data statistics
    data_dict_log = merge_fieldvals_and_main(data_dict_log,
                                             stats_dict_log,
                                             pepdep_xls_df)
    print(f"Counting proteins based on threshold {0.05} for raw data...")
    count_proteins(data_dict_vals,
                   comp_cols,
                   args[2],
                   meta_delim = args[4],
                   outpath = os.path.join(path, "raw"))
    print(f"Counting proteins based on threshold {0.05} for log data...")
    count_proteins(data_dict_log,
                   comp_cols,
                   args[2],
                   meta_delim = args[4],
                   outpath = os.path.join(path, "log"))
    # Then, use the dictionaries to make barplots.
    print(f"Plotting all bar plots to {args[2]}") 
    wrap_barplots(data_dict_vals,
                  cols_of_interest,
                  comp_cols,
                  outpath = os.path.join(path, "raw"))
    wrap_barplots(data_dict_log,
                  cols_of_interest,
                  comp_cols,
                  outpath = os.path.join(path, "log"))
    print("Complete!! :)\n")
    # and return None at the end.
    return None


#
#
############################################################################################################
#
#

if __name__ == "__main__":
    cmd_main()

#
#
############################################################################################################