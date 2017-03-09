#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to clean up and manage rsids.
"""
import numpy as np
import pandas as pd

NORM_FILE  = '/godot/dbsnp/rsid_lookup.pd'
DBSNP_FILE = '/godot/dbsnp/dbsnp149_lookup.snp_only.bed.rslookup.rs_sort.txt'


###############################################################################
#                            Clean Up rsid column                             #
###############################################################################


def clean_and_normalize_rsids(df, rsid_col='rsID', convert_df=NORM_FILE):
    """Clean and filter rsID column and normalize to current rsids only.

    After stripping whitespace, only rows with rsids that begin with 'rs' are
    kept (null rows are dropped) and then rsids are normalized to current
    rsids using normalize_rsids().

    Args:
        df (DataFrame):   The pandas dataframe to filter
        rsid_col (str):   The name os the rsid column
        convert_df (str): A file containing a normalization file from
                          make_norm_file()

    Returns:
        DataFrame: A matching dataframe with rsids corrected and bad rows
                   removed.
    """
    df = clean_rsids(df, rsid_col)
    df[rsid_col] = normalize_rsids(df[rsid_col], convert_df)
    return df


def clean_rsids(df, rsid_col='rsID'):
    """Strip whitespace from rsID column and return only valid rsIDs.

    After stripping whitespace, only rows with rsids that begin with 'rs' are
    kept (null rows are dropped).

    Args:
        df (DataFrame): The pandas dataframe to filter
        rsid_col (str): The name os the rsid column

    Returns:
        DataFrame: A matching dataframe with bad rsid rows removed.
    """
    orig = len(df)

    # Drop nan
    df  = df[df[rsid_col].notnull()].copy()
    nan = orig-len(df)

    # Stripe whitespace
    df[rsid_col] = df[rsid_col].str.strip()

    bad = (orig-nan)-len(df)

    df = df[df[rsid_col].str.startswith('rs')]

    print('Dropped {} bad rsids of {} total rows, {} for not starting with rs, {} for being nan'
          .format(nan+bad, orig, bad, nan))

    return df


def normalize_rsids(rsids, convert_df=NORM_FILE):
    """Update rsids to match current rsids from norm file.

    Args:
        rsids (Series):   A Series containing rsids to update.
        convert_df (str): A file containing a normalization file from
                          make_norm_file()
    Returns:
        Series: A series updated with latest rsids
    """
    assert isinstance(rsids, pd.core.series.Series)
    lookup = pd.read_pickle(convert_df)
    df     = pd.DataFrame(rsids)
    rcol   = df.columns[0]
    df     = pd.merge(df, lookup, how='left', left_on=rcol, right_index=True)
    return np.where(df.new_rs.notnull(), df.new_rs, df[rcol])


def make_norm_file(norm_file, outfile='rslookup.pd'):
    """Create a lookup dataframe from the RsMergeArch.bcp file from NCBI.

    Latest file here:
    ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz

    Args:
        norm_file (str): The above mentioned file path.
        outfile (str):   Where to write the pickled output to.
    """
    lookup = pd.read_csv(norm_file, sep='\t',
                         index_col=0, header=None)
    lookup.columns = ['new_rs', 'build_id', 'orien', 'create_time',
                      'last_update_time', 'rsCurrent', 'orien2Current',
                      'something']
    lookup.index.name = None

    lookup['old_rs'] = 'rs' + lookup.index.to_series().astype(str)
    lookup['new_rs'] = 'rs' + lookup.rsCurrent.astype(str)

    lookup.drop_duplicates(inplace=True)

    lookup.set_index('old_rs', drop=True, inplace=True)
    lookup = lookup[['new_rs']]

    lookup.to_pickle(outfile)
