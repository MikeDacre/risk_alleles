# -*- coding: utf-8 -*-
"""
SQLAlchemy database mappings for the risk allele database.
"""
import os as _os
from shutil import copy as _copy

from time import sleep as _sleep

from datetime import datetime as _dt

from urllib.request import urlopen as _req
from xml.dom.minidom import parseString as _ps

import pandas as _pd

from sqlalchemy import create_engine as _create_engine
from sqlalchemy import Column as _Column
from sqlalchemy import or_ as _or
from sqlalchemy import and_ as _and
from sqlalchemy import Date as _Date
from sqlalchemy import Float as _Float
from sqlalchemy import Index as _Index
from sqlalchemy import String as _String
from sqlalchemy import Integer as _Integer
from sqlalchemy import ForeignKey as _ForeignKey
from sqlalchemy import UniqueConstraint as _Unique

from sqlalchemy.orm import relationship as _relationship
from sqlalchemy.orm import sessionmaker as _sessionmaker
from sqlalchemy.ext.declarative import declarative_base as _base

Base = _base()

location = '/godot/risk_alleles/risk_alleles.db'

P_CUTOFF = 0.005


###############################################################################
#                              Table Declaration                              #
###############################################################################


class T(object):

    """A wrapper for all SQLAlchemy Tables."""

    class Trait(Base):

        """This class holds a lookup table of traits."""

        __tablename__ = 'traits'

        id      = _Column(_Integer, primary_key=True)
        trait   = _Column(_String, index=True, nullable=False)

        def __repr__(self):
            """Display summary stats."""
            return "Trait<id={id}, trait={trait}>".format(
                id=self.id, trait=self.trait
            )

        def __len__(self):
            """Number of SNPs."""
            return len(self.snps)

    class PMID(Base):

        """This class holds a lookup table of rsids."""

        __tablename__ = 'pmids'

        pmid    = _Column(_String, primary_key=True)
        title   = _Column(_String, nullable=True)
        pubdate = _Column(_Date, nullable=True)

        def __init__(self, pmid, add_study_data=True):
            """Get info from ID and set."""
            self.pmid = pmid
            if add_study_data:
                try:
                    title, pubdate = get_pubmed_info(pmid)
                    self.title   = title
                    self.pubdate = pubdate
                except:
                    pass

        def __repr__(self):
            """Display summary stats."""
            return "PMID<pmid={pmid}, title={title}, pubdate={pubdate}>".format(
                pmid=self.pmid, title=self.title, pubdate=self.pubdate
            )

        def __len__(self):
            """Number of SNPs."""
            return len(self.snps)

    class RiskAllele(Base):

        """This class holds the SQLAlchemy definition of the risk allele table.

        Always Populated Columns:
            id (int):          A unique ID
            pmid (int):        The pubmed ID of the study
            trait (str):       A string describing the trait of interest
            rsID (str):        The rsid of the SNP
            population (str):  The population the data is for
            risk_allele (str): The risk allele

        Optional Columns:
            chrom (str):      Chromosome string (hg19)
            position (int):   Position on chromosome (hg19)
            A1 (str):         The A1 allele
            A2 (str):         The A2 allele
            OR (float):       The odds ratio for A1
            B (float):        The beta (effect size) for A1
            OR_B (float):     Either the OR or B
            P (float):        p-value of the SNP-trait pair
            N (int):          Number of individuals supporting the SNP
            N_cases (int):    Number of of individuals in the case group
            N_controls (int): Number of of individuals in the control group
            source (str):     Which database this record is from
            grasp_id (int):   The ID of this SNP-trait pair in the GRASP database
            population_info (str): Extra population information
        """
        __tablename__ = 'risk_alleles'

        id              = _Column(_Integer, primary_key=True)
        pmid            = _Column(_String, _ForeignKey('pmids.pmid'), nullable=False)
        paper           = _relationship('PMID', backref='snps')
        trait_id        = _Column(_Integer, _ForeignKey('traits.id'), nullable=False)
        trait           = _relationship('Trait', backref='snps')
        rsID            = _Column(_String, index=True, nullable=False)
        population      = _Column(_String, index=True, nullable=False)
        risk_allele     = _Column(_String(1), nullable=False)
        chrom           = _Column(_String, index=True)
        position        = _Column(_Integer, index=True)
        A1              = _Column(_String)
        A2              = _Column(_String)
        OR              = _Column(_Float, index=True)
        B               = _Column(_Float)
        OR_B            = _Column(_Float)
        P               = _Column(_Float(precision=128), index=True)
        cases_MAF       = _Column(_Float)
        controls_MAF    = _Column(_Float)
        N               = _Column(_Integer)
        N_cases         = _Column(_Integer)
        N_controls      = _Column(_Integer)
        source          = _Column(_String)
        grasp_id        = _Column(_Integer)
        population_info = _Column(_String)

        pmid_rsid = _Index('pmid_rsid', pmid, rsID)
        loc       = _Index('loc', chrom, position)

        # Make pmid, trait, rsID, population, risk_allele unique
        _Unique(pmid, trait_id, rsID, population, risk_allele)

        def __repr__(self):
            """Better display of data."""
            out = ("RiskAllele<id={id}, pmid={pmid}, rsID={rsid}, trait={trait}, "
                   "risk_allele={risk_allele}>")
            return out.format(
                id=self.id, pmid=self.pmid, rsid=self.rsID, trait=self.trait.trait,
                risk_allele=self.risk_allele
            )

    class SigRiskAllele(Base):

        """This class holds the SQLAlchemy definition for significant risk allele.

        P-Value Cutoff = 0.005

        Always Populated Columns:
            id (int):          A unique ID
            pmid (int):        The pubmed ID of the study
            trait (str):       A string describing the trait of interest
            rsID (str):        The rsid of the SNP
            population (str):  The population the data is for
            risk_allele (str): The risk allele

        Optional Columns:
            chrom (str):      Chromosome string (hg19)
            position (int):   Position on chromosome (hg19)
            A1 (str):         The A1 allele
            A2 (str):         The A2 allele
            OR (float):       The odds ratio for A1
            B (float):        The beta (effect size) for A1
            OR_B (float):     Either the OR or B
            P (float):        p-value of the SNP-trait pair
            N (int):          Number of individuals supporting the SNP
            N_cases (int):    Number of of individuals in the case group
            N_controls (int): Number of of individuals in the control group
            source (str):     Which database this record is from
            grasp_id (int):   The ID of this SNP-trait pair in the GRASP database
            population_info (str): Extra population information
        """
        __tablename__ = 'sig_risk_alleles'

        id              = _Column(_Integer, primary_key=True)
        pmid            = _Column(_String, _ForeignKey('pmids.pmid'), nullable=False)
        paper           = _relationship('PMID', backref='sigsnps')
        trait_id        = _Column(_Integer, _ForeignKey('traits.id'), nullable=False)
        trait           = _relationship('Trait', backref='sigsnps')
        rsID            = _Column(_String, index=True, nullable=False)
        population      = _Column(_String, index=True, nullable=False)
        risk_allele     = _Column(_String(1), nullable=False)
        chrom           = _Column(_String, index=True)
        position        = _Column(_Integer, index=True)
        A1              = _Column(_String)
        A2              = _Column(_String)
        OR              = _Column(_Float, index=True)
        B               = _Column(_Float)
        OR_B            = _Column(_Float)
        P               = _Column(_Float(precision=128), index=True)
        cases_MAF       = _Column(_Float)
        controls_MAF    = _Column(_Float)
        N               = _Column(_Integer)
        N_cases         = _Column(_Integer)
        N_controls      = _Column(_Integer)
        source          = _Column(_String)
        grasp_id        = _Column(_Integer)
        population_info = _Column(_String)

        pmid_rsid = _Index('sig_pmid_rsid', pmid, rsID)
        loc       = _Index('sig_loc', chrom, position)

        # Make pmid, trait, rsID, population, risk_allele unique
        _Unique(pmid, trait_id, rsID, population, risk_allele)

        def __repr__(self):
            """Better display of data."""
            out = ("SigRiskAllele<id={id}, pmid={pmid}, rsID={rsid}, trait={trait}, "
                   "risk_allele={risk_allele}>")
            return out.format(
                id=self.id, pmid=self.pmid, rsid=self.rsID, trait=self.trait.trait,
                risk_allele=self.risk_allele
            )


###############################################################################
#                               Database Class                                #
###############################################################################


class DB(object):

    """A wrapper for the entire risk allele database."""

    location    = '/godot/risk_alleles/risk_alleles.db'

    basic_columns = [
        T.RiskAllele.id,
        T.RiskAllele.rsID,
        T.RiskAllele.population,
        T.Trait.trait,
        T.RiskAllele.risk_allele,
        T.RiskAllele.P,
        T.RiskAllele.pmid,
        T.PMID.title,
    ]

    extra_columns = [
        T.RiskAllele.chrom,
        T.RiskAllele.position,
        T.RiskAllele.A1,
        T.RiskAllele.A2,
        T.RiskAllele.OR,
        T.RiskAllele.B,
        T.RiskAllele.OR_B,
        T.RiskAllele.N,
        T.RiskAllele.N_cases,
        T.RiskAllele.N_controls,
        T.RiskAllele.source,
        T.RiskAllele.grasp_id,
        T.RiskAllele.population_info,
    ]
    basic_columns_sig = [
        T.SigRiskAllele.id,
        T.SigRiskAllele.rsID,
        T.SigRiskAllele.population,
        T.Trait.trait,
        T.SigRiskAllele.risk_allele,
        T.SigRiskAllele.P,
        T.SigRiskAllele.pmid,
        T.PMID.title,
    ]
    extra_columns_sig = [
        T.SigRiskAllele.chrom,
        T.SigRiskAllele.position,
        T.SigRiskAllele.A1,
        T.SigRiskAllele.A2,
        T.SigRiskAllele.OR,
        T.SigRiskAllele.B,
        T.SigRiskAllele.OR_B,
        T.SigRiskAllele.N,
        T.SigRiskAllele.N_cases,
        T.SigRiskAllele.N_controls,
        T.SigRiskAllele.source,
        T.SigRiskAllele.grasp_id,
        T.SigRiskAllele.population_info,
    ]

    def __init__(self, loc=None, debug=False):
        """Attach to a database, create if does not exist."""
        if not loc:
            loc = self.location
        self.location = _os.path.abspath(location)
        kwargs = {'echo': True} if debug else {'echo': False}
        self.engine   = _create_engine('sqlite:///{}'.format(self.location),
                                       **kwargs)
        if not _os.path.isfile(self.location):
            self.create_database()
        self.columns = list(T.RiskAllele.__table__.columns.keys())

    ##########################################################################
    #                           Basic Connectivity                           #
    ##########################################################################

    def get_session(self):
        """Return engine, session for this database."""
        Session = _sessionmaker(bind=self.engine)
        return Session()

    @property
    def session(self):
        """Simple wrapper for a new session."""
        return self.get_session()

    ##########################################################################
    #                              Adding Data                               #
    ##########################################################################

    def add_new_study(self, df, confirm=True):
        """Take a dataframe with the same columns as us and write it to self.

        Note: Will append all data to the database, which means it will not
              clobber any existing data, so you can *easily* end up with
              duplicated data.

        Also coerces some common incorrect columns to the correct format and
        drops non-allowed columns
        """
        # Check if study in DB already
        your_pmids = df.pmid.unique()
        our_pmids  = self.pmids
        for pmid in your_pmids:
            if pmid in our_pmids and confirm:
                ans = input('PMID {} already in DB, add anyway? [y/N] '
                            .format(pmid))
                if not ans.upper().startswith('Y'):
                    print('Aborting')
                    return

        df = df.copy()

        traits = [{'trait': t} for t in df.trait.unique() if t not in self.traits]
        if traits:
            print('Adding missing traits')
            conn = self.engine.connect()
            ins = T.Trait.__table__.insert()
            conn.execute(ins, traits)
            conn.close()
        traits = _pd.DataFrame.from_dict(self.trait_ids, orient='index')
        traits.columns = ['trait_id']

        print('Adding traits to dataframe')
        df = _pd.merge(df, traits, how='left', left_on='trait',
                       right_index=True)
        df = df.drop('trait', axis=1)

        print('Adding studies')
        for pmid in [str(i) for i in df.pmid.unique()]:
            if pmid not in self.pmids:
                print(pmid)
                self.add_pmid(pmid)
                _sleep(0.2)

        # Correct columns
        print('Cleaning columns')
        lookup = {
            'OR/B':         'OR_B',
            'start':        'position',
            'cases_maf':    'cases_MAF',
            'control_maf':  'controls_MAF',
            'controls_maf': 'controls_MAF',
            'N_CASES':      'N_cases',
            'N_CONTROLS':   'N_controls',
        }
        new_cols = []
        for col in df.columns:
            new_cols.append(lookup[col] if col in lookup else col)
        df.columns = new_cols

        # Find non-matching columns
        ok_cols  = []
        bad_cols = []
        for col in df.columns:
            if col in self.columns:
                ok_cols.append(col)
            else:
                bad_cols.append(col)
        if bad_cols and confirm:
            ans = input('Columns {} are not allowed in the db and will be '
                        .format(bad_cols) +
                        'dropped. Continue writing to the DB? [y/N] ')
            if not ans.upper().startswith('Y'):
                print('Aborting')
                return

        # Filter columns
        df = df[ok_cols]

        # Add to db using pandas syntax
        print('Adding data')
        df.to_sql('risk_alleles',
                  self.engine,
                  chunksize=100000,
                  index=False,
                  if_exists='append')

        # Add significant SNPs to sig table
        df2 = df[df.P < P_CUTOFF].copy()
        print('Adding signficiant SNP table')
        df2.to_sql('sig_risk_alleles',
                  self.engine,
                  chunksize=100000,
                  index=False,
                  if_exists='append')

    def add_pmid(self, pmid):
        """Check if pmid is in DB, if not, add. Return ID."""
        if pmid not in self.pmids:
            print('Adding pubmed id {}'.format(pmid))
            session = self.get_session()
            new_pmid = T.PMID(pmid=pmid, add_study_data=True)
            session.add(new_pmid)
            session.commit()
        return pmid

    def add_trait(self, trait):
        """Check if trait is in DB, if not, add. Return ID."""
        session = self.get_session()
        if trait not in self.traits:
            print('Adding trait {}'.format(trait))
            new_trait = T.Trait(trait=trait)
            session.add(new_trait)
            session.commit()
        i = session.query(T.Trait.id).filter(T.Trait.trait == trait).first()
        return i[0]

    ##########################################################################
    #                                Querying                                #
    ##########################################################################

    def query(self, *args, sig_only=False, **kwargs):
        """Wrapper for the SQLAlchemy query method of session.

        If there are no args, then the whole table is returned, sig_only
        modifies if the significant allele table or the whole table.

        Args:
            sig_only (bool): Limit to significant alleles table only.
            *args, **kwargs: Any arguments allowed by session.query.
        """
        if not args:
            rtable = T.SigRiskAllele if sig_only else T.RiskAllele
            args = (rtable,)
        session = self.get_session()
        return session.query(*args, **kwargs)

    def get_dataframe(self, pval=None, rsid=None, trait=None, population=None,
                      position=None, limit=None, sig_only=False, extra=False):
        """Filter database and return results as a dataframe.

        All arguments are optional and additive, no arguments will return the
        whole database.

        Note that using multiple arguments at the same time can make the
        operation very, very slow, try to use as few as possible.

        Args:
            pval (float):          Filter by pvalue less than this number.
            rsid (str/list):       Filter by rsid, can be a list of rsids.
            trait (str/list):      Filter by trait name, can be a list.
            population (str/list): Filter by population, can be a list.
            position (str, list):  Filter by location, syntax is important,
                                   must be:
                                       chr:start-end
                                   End is optional, both numbers are base-0.
                                   You can provide a list of these, but each
                                   must match the syntax requirements.
            limit (int):           Number of results to limit to.
            sig_only (bool):       Only return results in the significant table.
            extra (bool):          Add additional information columns to output.

        Returns:
            DataFrame: A pandas dataframe of the results.
        """
        # Choose columns
        rtable = T.SigRiskAllele if sig_only else T.RiskAllele
        columns = self.basic_columns_sig if sig_only else self.basic_columns
        if extra:
            columns += self.extra_columns_sig if sig_only else self.extra_columns

        # Get the initial query object
        query = self.query(*columns)

        # Run filters
        if pval:
            assert isinstance(pval, (float, int))
            query = query.filter(rtable.P < pval)
        if rsid:
            if isinstance(rsid, str):
                query = query.filter(rtable.rsID == rsid)
            else:
                query = query.filter(rtable.rsID.in_(listify(rsid)))
        if trait:
            if isinstance(trait, int):
                query = query.filter(rtable.trait_id == trait)
            elif isinstance(trait, str):
                query = query.filter(T.Trait.trait == trait)
                query = query.filter(rtable.trait_id == T.Trait.id)
            else:
                trait = listify(trait)
                if isinstance(trait, int):
                    query = query.filter(rtable.trait.in_(trait))
                else:
                    query = query.filter(T.Trait.trait.in_(trait))
                    query = query.filter(rtable.trait_id.in_(T.Trait.id))
        if population:
            if isinstance(population, str):
                query = query.filter(rtable.population == population)
            else:
                query = query.filter(rtable.population.in_(
                    listify(population)
                ))
        if position:
            if isinstance(position, str):
                chrom, start, end = parse_position(position)
                query = query.filter(rtable.chrom == chrom)
                if start:
                    if end:
                        query = query.filter(
                            rtable.position.between(start, end)
                        )
                    else:
                        query = query.filter(rtable.position == start)
            else:
                positions = listify(position)
                or_list = []
                for pos in positions:
                    chrom, start, end = parse_position(pos)
                    if start:
                        if end:
                            or_list.append(
                                _and(
                                    rtable.chrom == chrom,
                                    rtable.position.between(start, end)
                                )
                            )
                        else:
                            or_list.append(
                                _and(
                                    rtable.chrom == chrom,
                                    rtable.position == start
                                )
                            )
                    else:
                        or_list.append(rtable.chrom == chrom)
                query = query.filter(_or(*or_list))

        # Set limits
        if limit:
            query = query.limit(int(limit))

        # Get DataFrame
        df = _pd.read_sql_query(query.statement, self.engine, index_col='id')
        df.index.name = None
        return df

    @property
    def significant(self):
        """Return all significant SNPs as a DataFrame."""
        session = self.get_session()
        return _pd.read_sql_query(
            session.query(*self.basic_columns_sig).filter(
                T.SigRiskAllele.trait_id == T.Trait.id
            ).filter(
                T.SigRiskAllele.pmid == T.PMID.pmid
            ).statement,
            self.engine,
            index_col='id'
        )

    @property
    def pmids(self):
        """Return a list of pubmed ids in self."""
        return [i[0] for i in self.query(T.PMID.pmid).all()]

    @property
    def traits(self):
        """Return a list of traits."""
        return sorted([
            t[0] for t in self.query(
                T.Trait.trait,
            ).all()
        ])

    @property
    def trait_ids(self):
        """Return a dictionary of trait=>id."""
        return {
            t: i for t, i in self.query(
                T.Trait.trait,
                T.Trait.id
            ).all()
        }

    ##########################################################################
    #                          Maintenance Methods                           #
    ##########################################################################

    def create_database(self, confirm=True):
        """Create the db from scratch.

        Note: Will DELETE the current database!!!
        """
        if confirm:
            ans = input('Are you sure you want to erase and recreate the db? ')
            if not ans.upper().startswith('Y'):
                print('Aborting')
                return
        if _os.path.exists(self.location):
            _os.remove(self.location)
        Base.metadata.create_all(self.engine)
        print('Done')

    def backup(self):
        """Backup the database."""
        timestamp = '{:%Y-%m-%d %H:%M:%S}'.format(_dt.now())
        new = '{}.{}'.format(self.location, timestamp)
        _copy(self.location, new)
        print('Backup at {}'.format(new))

    ##########################################################################
    #                               Internals                                #
    ##########################################################################

    def __getitem__(self, x):
        """Quick access to rsids."""
        if isinstance(x, str):
            return self.query().filter(T.RiskAllele.rsID == x).all()
        try:
            x = list(x)
        except TypeError:
            raise TypeError('Cannot lookup {}'.format(x))
        return self.query().filter(T.RiskAllele.rsID.in_(x)).all()

    def __len__(self):
        """Print the length."""
        session = self.get_session()
        return session.query(T.RiskAllele).count()

    def __repr__(self):
        """Basic information about self."""
        return 'RiskAlleles<{location}>'.format(location=self.location)

    def __str__(self):
        """Basic information about self."""
        return 'RiskAlleles<{location}>'.format(location=self.location)



###############################################################################
#                              Helper Functions                               #
###############################################################################


def get_pubmed_info(pmid):
    """Get title and pubdate from pubmed id.

    Returns:
        (str, datetime): title, publication date
    """
    print('Getting study data for {}'.format(pmid))

    data = _ps(
        _req(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
            '?db=pubmed&retmode=xml&id={id}'.format(id=pmid),
            timeout=60
        ).read()
    )

    article = data.getElementsByTagName(
        'PubmedArticle'
    )[0].getElementsByTagName(
        'MedlineCitation'
    )[0].getElementsByTagName(
        'Article'
    )[0]

    title = article.getElementsByTagName(
        'ArticleTitle'
    )[0].childNodes[0].nodeValue

    date = data.getElementsByTagName(
        'Journal'
    )[0].getElementsByTagName(
        'JournalIssue'
    )[0].getElementsByTagName(
        'PubDate'
    )[0]

    year  = date.getElementsByTagName('Year')[0].childNodes[0].nodeValue
    month = date.getElementsByTagName('Month')[0].childNodes[0].nodeValue
    day   = date.getElementsByTagName('Day')[0].childNodes[0].nodeValue
    pubdate = _dt.strptime(
        '{year} {month} {day}'.format(year=year, month=month, day=day),
        '%Y %b %d'
    )
    return title, pubdate


def listify(x):
    """Turn an iterable into a list if possible."""
    if isinstance(x, (str, int, float)):
        return [x]
    return list(x)


def parse_position(x):
    """Parse chrom:start-end into chrom, start, end."""
    if ':' in x:
        chrom, rng = x.split(':')
    else:
        return x, None, None
    if '-' in chrom:
        start, end = rng.split('-')
        return chrom, int(start), int(end)
    return chrom, int(rng), None
