#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A bunch of little functions to work with dbSNP files quickly.

Simple workflow:
    1. Download dbSNP from NCBI
    2. Concatenate all chromosome bed files into a single bed (eg. bed149.bed)
    3. Run make_lookup_tables(): dbsnp create -o dbsnp149 <dbsnp149.bed>
    4. While that is running, dump all your rsids to a file
    5. Run join_rsid(): dbsnp rsids.txt dbsnp149.rslookup.rs_sort.txt rsids.joined.txt

The same thing can be done for locations, just dump your locations to a file
with the locations filtered as chr.position. Then provide either
dbsnp149.rslookup.start_sort.txt or dbsnp149.rslookup.end_sort.txt as the
dnsnp_file, use the 'end' file if your data is base-1 coded.

This file can also be imported and the functions used directly, this provides
the added advantage of being able to pass Series and get back DataFrames
instead of having to work with the files.

Note: This script will only work on linux systems, it makes use of the linux
tools sort, awk, cat, zcat, and join; these tools work differently on different
systems. Be sure to check all outputs for sanity.

Also, dbSNP is huge, and these commands are run in serial, so they will take
multiple minutes or more to complete, particularly when creating the lookup
files.
"""
import os as _os
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse
from random import randint as _rand
from subprocess import check_call as _call
from subprocess import STDOUT as _STDOUT
from subprocess import CalledProcessError as _call_err

import pandas as pd

__all__ = ['make_lookup_tables', 'join_rsid', 'join_location']


###############################################################################
#                        dbSNP Manipulation Functions                         #
###############################################################################


def join_rsid(rsids, dbsnp_file, outfile, sort=True, as_df=False):
    """Use linux join to create a lookup table of rsids.

    Args:
        rsids (str/list): List of rsids as a file name (string), list of rsid,
                          or Series.
        dbsnp_file (str): The dbsnp lookup file from make_lookup_tables.
                          should be the .rslookup.rs_sort.txt file (zipped ok)
        outfile (str):    Name of outfile to write to
        sort (bool):      Pre sort the rsids
        as_df (bool):     Return a dataframe

    Writes:
        A tab separated table of rsid, chrom, start, end for all rsids.

    Returns:
        DataFrame: Dataframe of written table, only returned if as_df is true.
    """
    if isinstance(rsids, pd.core.series.Series):
        rsids = rsids.tolist()

    if isinstance(rsids, (list, tuple, set)):
        rsids = sorted(list(set(rsids)))
        tmpfile = outfile + '.rsids.tmp'
        with open(tmpfile, 'w') as fout:
            fout.write('\n'.join(rsids))
        rsids = tmpfile
    else:
        tmpfile = None

    rsids   = _os.path.abspath(rsids)
    outfile = _os.path.abspath(outfile)

    if sort:
        print('Sorting')
        cat = 'zcat' if rsids.endswith('gz') else 'cat'
        tmpfile = 'tmpsort_{}'.format(_rand(1000,20000))
        script = r"""{cat} {rsids} | sort > {tmp}; mv {tmp} {rsids}"""
        _call(script.format(cat=cat, rsids=rsids, tmp=tmpfile), shell=True)

    print('Joining')
    script = r"""join {rsids} {dbsnp} > {outfile}"""
    try:
        _call(
            script.format(rsids=rsids, dbsnp=dbsnp_file, outfile=outfile),
            stderr=_STDOUT, shell=True, universal_newlines=True
        )
    except _call_err as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        raise exc
    print('Done, file {} has the joined list'.format(outfile))

    if as_df:
        print('Getting DataFrame')
        try:
            df = pd.read_csv(outfile, sep=' ', header=None, index_col=0)
        except pd.io.common.EmptyDataError:
            print('Joined file empty, skipping')
            return None
        df.index.name = None
        df.columns = ['chrom', 'start', 'end']
        return df


def join_location(locs, dbsnp_file, outfile, sort=True, as_df=False):
    """Use linux join to create a lookup table of rsids.

    Args:
        locs (str/list):  File name of a file of chr.position, of a list of the
                          same. Can also give a Series of the same, or a
                          DataFrame with the columns 'chrom' and 'position'.
        dbsnp_file (str): The dbsnp lookup file from make_lookup_tables.
                          should be the .rslookup.rs_sort.txt file (zipped ok)
        outfile (str):    Name of outfile to write to
        sort (bool):      Pre sort the rsids
        as_df (bool):     Return a dataframe

    Writes:
        A tab separated table of chrom.pos, rsid for all locations.

    Returns:
        DataFrame: Dataframe of written table, only returned if as_df is true.
    """
    if isinstance(locs, pd.core.series.Series):
        locs = sorted(list(set(locs.tolist())))

    if isinstance(locs, pd.core.frame.DataFrame):
        locs = locs.chrom.astype(str) + '.' + locs.position.astype(str)
        locs = sorted(list(set(locs.tolist())))

    if isinstance(locs, (list, tuple, set)):
        tmpfile = outfile + 'locs.tmp'
        with open(tmpfile, 'w') as fout:
            fout.write('\n'.join(locs))
        locs = tmpfile
    else:
        tmpfile = None

    locs    = _os.path.abspath(locs)
    outfile = _os.path.abspath(outfile)

    if sort:
        print('Sorting')
        cat = 'zcat' if locs.endswith('gz') else 'cat'
        script = r"""{cat} {locs} | sort -k1,1 > tmp45876; mv tmp45876 {locs}"""
        _call(script.format(cat=cat, locs=locs), shell=True)

    print('Joining')
    script = r"""join {locs} {dbsnp} > {outfile}"""
    _call(script.format(locs=locs, dbsnp=dbsnp_file, outfile=outfile),
          shell=True)

    if tmpfile:
        _os.remove(tmpfile)

    if as_df:
        print('Getting DataFrame')
        try:
            lookup = pd.read_csv(outfile, sep=' ', header=None, index_col=0)
        except pd.io.common.EmptyDataError:
            print('Joined file empty, skipping')
            return None
        lookup.index.name = None
        lookup['chrom'], lookup['position'] = lookup.index.to_series().str.split('.', 1).str
        lookup.columns = ['rsid', 'chrom', 'position']
        lookup = lookup[['chrom', 'position', 'rsid']]
        return lookup


def make_lookup_tables(dbsnp_file, outname='dbsnp_filtered',
                       skip_filter=False):
    """Create a series of easy lookup tables from dbSNP files."""
    filtered_bed = '{}.snp_only.bed'.format(outname)
    if not skip_filter:
        print('Making filtered snp_only bed file')
        with open_zipped(dbsnp_file) as fin, open(filtered_bed, 'w') as fout:
            for line in fin:
                if line.startswith('track'):
                    continue
                chrom, start, end, name, _, strand = line.rstrip().split('\t')
                if int(end)-int(start) > 1:
                    continue
                fout.write('\t'.join([chrom, start, end, name, strand]) + '\n')
    # Create a lookup file sorted by rsid
    print('Making rsid lookup table')
    script = r"""cat {name} | awk '{{print $4 "\t" $1 "\t" $2 "\t" $3}}' | sort -k1,1 > {name}.rslookup.rs_sort.txt"""
    _call(script.format(name=filtered_bed), shell=True)
    # Create a location lookup file sorted by the start location (base-1)
    print('Making start location lookup table')
    script = r"""cat {name} | awk '{{print $1"."$2 "\t" $4}}' | sort -k1,1 > {name}.rslookup.start_sort.txt"""
    _call(script.format(name=filtered_bed), shell=True)
    # Create a location lookup file sorted by the end location (base-1)
    print('Making end location lookup table')
    script = r"""cat {name} | awk '{{print $1"."$3 "\t" $4}}' | sort -k1,1 > {name}.rslookup.end_sort.txt"""
    _call(script.format(name=filtered_bed), shell=True)
    print('Done')


###############################################################################
#                              Helper Functions                               #
###############################################################################


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of zipped or not.

    Text mode enforced for compatibility with python2.
    """
    mode   = mode[0] + 't'
    p2mode = mode
    if hasattr(infile, 'write'):
        return infile
    if isinstance(infile, str):
        if infile.endswith('.gz'):
            return _gzip.open(infile, mode)
        if infile.endswith('.bz2'):
            if hasattr(_bz2, 'open'):
                return _bz2.open(infile, mode)
            else:
                return _bz2.BZ2File(infile, p2mode)
        return open(infile, p2mode)


def make_arg_parser():
    """Parse command line arguments."""
    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(title='mode', dest='mode')

    # Create lookup tables
    create = subparsers.add_parser('create', aliases=['c'],
                                   help="Create the lookup tables")
    create.add_argument('dbsnp_file', help="dbSNP file, raw download.")
    create.add_argument('-o', '--outname',
                        help="A name prefix for the output files.")

    # rsid lookup
    rsids = subparsers.add_parser('rsids', aliases=['rsid', 'r'],
                                  help="Get locations by rsid")
    rsids.add_argument('rsid_file', help='File of rsids to join')
    rsids.add_argument('dbsnp_file', help="dbSNP rsid lookup file.")
    rsids.add_argument('outfile', help="File to write output to.")
    rsids.add_argument('--no-sort', action='store_false',
                       help="rsid_file is already sorted")

    # location lookup
    locs = subparsers.add_parser('location', aliases=['locs', 'loc', 'l'],
                                 help="Get rsids by location")
    locs.add_argument('loc_file', help="File of chr.location")
    locs.add_argument('dbsnp_file',
                      help="dbSNP location lookup file be sure to chose " +
                      "the end location file if your data is base-1")
    locs.add_argument('outfile', help="File to write output to.")
    locs.add_argument('--no-sort', action='store_false',
                      help="loc_file is already sorted")

    return parser


###############################################################################
#                              Running Directly                               #
###############################################################################


def main(argv=None):
    """Run as a script."""
    # Command line args
    parser = make_arg_parser()
    if not argv:
        argv = _sys.argv[1:]
    args = parser.parse_args(argv)

    if args.mode == 'create':
        make_lookup_tables(args.dbsnp_file, args.outname)
    if args.mode == 'rsids':
        join_rsid(args.rsid_file, args.dbsnp_file, args.outfile, args.no_sort)
    if args.mode == 'location':
        join_location(args.loc_file, args.dbsnp_file, args.outfile, args.no_sort)


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
