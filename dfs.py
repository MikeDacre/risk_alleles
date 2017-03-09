#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Super functions to handle doing several things to a dataframe at once.
"""
import os as _os
from tempfile import NamedTemporaryFile as _tf

import pandas as pd

import dbsnp
import rsids

NORM_FILE  = rsids.NORM_FILE
DBSNP_FILE = rsids.DBSNP_FILE

###############################################################################
#                     Super Wrapper for Study DataFrames                      #
###############################################################################


def run_all(df, get_locs=False, rsid_col='rsID', risk_allele_col='risk_allele',
            columns=None, convert_df=NORM_FILE, dbsnp_file=DBSNP_FILE):
    """Fully clean and normalize df and add chromosome locations.

    Runs all clean functions to remove excess whitespace in every column,
    coerce risk alleles into A, T, G, or C, and drop all others.

    Converts rsIDs to most recent version and drop ones without an initial
    'rs'.

    Drops all columns in columns, defaults to:
        ['pmid', 'trait', 'rsID', 'population', 'risk_allele']

    Attempts to add locations for all known rsIDs.

    Returns a similar dataframe, cleaned up and with locations added.

    Args:
        df (DataFrame):           The pandas dataframe to filter
        get_locs (bool):          Get locations after cleaning
        rsid_col (str):           The name os the rsid column
        risk_allele_column (str): The name of the risk allele column.
        columns (list):           A list of columns names to drop if null.
        convert_df (str):         A file containing a normalization file from
                                  make_norm_file()
        dbsnp_file (str):         The dbsnp lookup file from make_lookup_tables
                                  should be the .rslookup.rs_sort.txt file

    Returns:
        DataFrame: A matching dataframe with rsids corrected, bad rows
                   (malformatted rsids) removed, and chrom, start, end columns
                   added.
    """
    print('Dropping null columns')
    df = drop_null(df, columns)
    df = df.copy()  # Avoid annoying warnings about slicing
    if get_locs:
        df = clean_df_add_locations(df, rsid_col, risk_allele_col, convert_df,
                                    dbsnp_file)
    else:
        df = clean_df(df, rsid_col, risk_allele_col, convert_df)

    return df


def clean_df(df, rsid_col='rsID', risk_allele_col='risk_allele',
             convert_df=NORM_FILE):
    """Clean a dataframe of rsids and risk_alleles.

    Args:
        df (DataFrame):           The pandas dataframe to filter
        rsid_col (str):           The name os the rsid column
        risk_allele_column (str): The name of the risk allele column.
        convert_df (str):         A file containing a normalization file from
                                  make_norm_file()

    Returns:
        DataFrame: A matching DataFrame with rsids and risk alleles corrected.
    """
    print('Stripping whitespace')
    df = clean_df_whitespace(df)
    print('Cleaning risk alleles')
    df = clean_risk_alleles(df, risk_allele_col)
    print('Normalizing rsIDs')
    df = rsids.clean_and_normalize_rsids(df, rsid_col, convert_df)
    print('Force PMIDs to string')
    df['pmid'] = df.pmid.astype(str)
    return df


def clean_df_add_locations(df, rsid_col='rsID', risk_allele_col='risk_allele',
                           convert_df=NORM_FILE, dbsnp_file=DBSNP_FILE):
    """Fully clean and normalize df and add chromosome locations.

    Runs clean_and_normalize_rsids() and then fetches location info.

    Args:
        df (DataFrame):           The pandas dataframe to filter
        rsid_col (str):           The name os the rsid column
        risk_allele_column (str): The name of the risk allele column.
        convert_df (str):         A file containing a normalization file from
                                  make_norm_file()
        dbsnp_file (str):         The dbsnp lookup file from make_lookup_tables
                                  should be the .rslookup.rs_sort.txt file

    Returns:
        DataFrame: A matching dataframe with rsids corrected, bad rows
                   (malformatted rsids) removed, and chrom, start, end columns
                   added.
    """
    clean_df(df, rsid_col, risk_allele_col, convert_df)
    print('Getting locations')
    return clean_rsids_add_locations(df, rsid_col, convert_df, dbsnp_file)


###############################################################################
#                           Fetch Location Columns                            #
###############################################################################


def clean_rsids_add_locations(df, rsid_col='rsID',
                              convert_df=NORM_FILE, dbsnp_file=DBSNP_FILE):
    """Fully clean and normalize rsids and add chromosome location to df.

    Runs normalize_rsids() and then fetches location info.

    Args:
        df (DataFrame):   The pandas dataframe to filter
        rsid_col (str):   The name os the rsid column
        convert_df (str): A file containing a normalization file from
                          make_norm_file()
        dbsnp_file (str): The dbsnp lookup file from make_lookup_tables.
                          should be the .rslookup.rs_sort.txt file (zipped ok)

    Returns:
        DataFrame: A matching dataframe with rsids corrected, bad rows
                   (malformatted rsids) removed, and chrom, start, end columns
                   added.
    """
    df[rsid_col] = rsids.normalize_rsids(df[rsid_col], convert_df)
    tf = get_tempfile_name()
    locations = dbsnp.join_rsid(df[rsid_col], dbsnp_file, tf, as_df=True)
    _os.remove(tf)
    if locations is not None:
        print('{} of {} rows have new locations, merging'
              .format(len(locations), len(df)))
        return pd.merge(df, locations, how='left', left_on=rsid_col,
                        right_index=True)
    else:
        print('No rows have locations in dbSNP file')
        return df


###############################################################################
#                                  Clean DF                                   #
###############################################################################


def clean_df_whitespace(df):
    """Clean all columns of whitespace."""
    for column in df.columns:
        if df[column].dtype.name is 'object':
            df[column] = df[column].astype(str).strip()
    return df


def clean_risk_alleles(df, risk_allele_column='risk_allele'):
    """Coerce risk alleles into A,G,C,T or drop on failure.

    Args:
        df (DataFrame):           Any pandas DataFrame.
        risk_allele_column (str): The name of the risk allele column.

    Returns:
        DataFrame: Same DataFrame as df, with risk alleles cleaned and dropped.
    """
    df[risk_allele_column] = df[risk_allele_column].str.upper()
    df[risk_allele_column] = df[risk_allele_column].str.strip()
    orig = len(df)
    df = df[df[risk_allele_column].isin(list('ATGC'))]
    final = len(df)
    print('{} bad risk alleles dropped for not being A/T/G/C'
          .format(orig-final))
    return df


def drop_null(df, columns=None):
    """Drop all rows will null entries in any of columns."""
    if not columns:
        columns = ['pmid', 'trait', 'rsID', 'population', 'risk_allele']
    print('Null rows:')
    for i in columns:
        print('\t', i, len(df[df[i].isnull()]))
    orig = len(df)
    df.dropna(subset=columns, inplace=True)
    final = len(df)
    print('Dropped {} null rows'.format(orig-final))
    return df


###############################################################################
#                              Helper Functions                               #
###############################################################################


def get_tempfile_name():
    """Hack the tempfile system for a name."""
    t = _tf('w')
    tf = t.name
    t.close()
    return tf
