# -*- coding: utf-8 -*-
"""
SQLAlchemy database mappings for the risk allele database.
"""
import os as _os
import sys as _sys
from shutil import copy as _copy

from time import sleep as _sleep

from datetime import datetime as _dt

from urllib.request import urlopen as _req
from xml.dom.minidom import parseString as _ps

import pandas as _pd

from sqlalchemy import create_engine as _create_engine
from sqlalchemy import Column as _Column
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


###############################################################################
#                              Table Declaration                              #
###############################################################################


class RiskAllele(Base):

    """This class holds the SQLAlchemy definition of the risk allele database.

    Always Populated Columns:
        id (int):          A unique ID
        pmid (int):        The pubmed ID of the study
        trait (str):       A string describing the trait of interest
        rsID (str):        The rsid of the SNP
        population (str):  The population the data is for
        risk_allele (str): The risk allele

    Optional Columns:
        chrom (str):      Chromosome string (hg19)
        pos (int):        Position on chromosome (hg19)
        A1 (str):         The A1 allele
        A2 (str):         The A2 allele
        OR (float):       The odds ratio for A1
        B (float):        The beta (effect size) for A1
        OR_B (float):     Either the OR or B
        P (float):        p-value of the SNP-trait pair
        N (int):          Number of individuals supporting the SNP
        N_Cases (int):    Number of of individuals in the case group
        N_Controls (int): Number of of individuals in the control group
        source (str):     Which database this record is from
        grasp_id (int):   The ID of this SNP-trait pair in the GRASP database
        population_info (str): Extra population information
    """
    __tablename__ = 'risk_alleles'

    id              = _Column(_Integer, primary_key=True)
    pmid            = _Column(_String, _ForeignKey('pmids.pmid'), nullable=False)
    trait_id        = _Column(_Integer, _ForeignKey('traits.id'), nullable=False)
    trait           = _relationship('Trait', back_populates='snps')
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


class Trait(Base):

    """This class holds a lookup table of traits."""

    __tablename__ = 'traits'

    id    = _Column(_Integer, primary_key=True)
    trait = _Column(_String, index=True, nullable=False)
    snps  = _relationship('RiskAllele', back_populates='trait')

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
    snps    = _relationship('RiskAllele')

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


###############################################################################
#                               Database Class                                #
###############################################################################


class RiskAlleles(object):

    """A wrapper for the entire risk allele database."""

    snp_table   = RiskAllele
    trait_table = Trait
    pmid_table  = PMID
    columns     = list(snp_table.__table__.columns.keys())
    location    = '/godot/risk_alleles/risk_alleles.db'

    def __init__(self, loc=None):
        """Attach to a database, create if does not exist."""
        if not loc:
            loc = self.location
        self.location = _os.path.abspath(location)
        self.engine   = _create_engine('sqlite:///{}'.format(self.location))
        if not _os.path.isfile(self.location):
            self.create_database()

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
            ins = self.trait_table.__table__.insert()
            conn.execute(ins, traits)
        traits = _pd.DataFrame.from_dict(self.traits, orient='index')
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

    def add_pmid(self, pmid):
        """Check if pmid is in DB, if not, add. Return ID."""
        if pmid not in self.pmids:
            print('Adding pubmed id {}'.format(pmid))
            session = self.get_session()
            new_pmid = PMID(pmid=pmid, add_study_data=True)
            session.add(new_pmid)
            session.commit()
        return pmid

    def add_trait(self, trait):
        """Check if trait is in DB, if not, add. Return ID."""
        session = self.get_session()
        if trait not in self.traits:
            print('Adding trait {}'.format(trait))
            new_trait = Trait(trait=trait)
            session.add(new_trait)
            session.commit()
        i = session.query(Trait.id).filter(Trait.trait == trait).first()
        return i[0]

    ##########################################################################
    #                                Querying                                #
    ##########################################################################

    def query(self, *args):
        """Wrapper for the SQLAlchemy query method of session.

        Args:
            *args: Any arguments allowed by session.query. If blank, self.table
                   is used, which will return the whole database. To limit by
                   columns, simply pass columns: query(self, self.table.rsID)
                   would return only a list of rsids.
        """
        if not args:
            args = (self.snp_table,)
        session = self.get_session()
        return session.query(*args)

    @property
    def pmids(self):
        """Return a list of pubmed ids in self."""
        return [i[0] for i in self.query(self.pmid_table.pmid).all()]

    @property
    def traits(self):
        """Return a dictionary of trait=>id."""
        return {
            t: i for t, i in self.query(
                self.trait_table.trait,
                self.trait_table.id
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
            return self.query().filter(self.snp_table.rsID == x).all()
        try:
            x = list(x)
        except TypeError:
            raise TypeError('Cannot lookup {}'.format(x))
        return self.query().filter(self.snp_table.rsID.in_(x)).all()

    def __len__(self):
        """Print the length."""
        session = self.get_session()
        return session.query(self.snp_table).count()

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
