import numpy as np
import pandas as pd

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, Float, String, Index, Boolean
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

location = 'combined.db'


class RiskAlleles(Base):

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

    id              = Column(Integer, primary_key=True)
    pmid            = Column(Integer, index=True, nullable=False)
    trait           = Column(String, index=True, nullable=False)
    rsID            = Column(String, index=True, nullable=False)
    population      = Column(String, index=True, nullable=False)
    risk_allele     = Column(String(1), nullable=False)
    chrom           = Column(String, index=True)
    pos             = Column(Integer, index=True)
    A1              = Column(String)
    A2              = Column(String)
    OR              = Column(Float, index=True)
    B               = Column(Float)
    OR_B            = Column(Float)
    P               = Column(Float(precision=128), index=True)
    N               = Column(Integer)
    N_Cases         = Column(Integer)
    N_Controls      = Column(Integer)
    source          = Column(String)
    grasp_id        = Column(Integer)
    population_info = Column(String)

    pmid_rsid = Index('pmid_rsid', pmid, rsID)
    loc       = Index('loc', chrom, pos)


def get_session():
    """Return engine, session for this database."""
    engine  = create_engine('sqlite:///{}'.format('combined.db'))
    Session = sessionmaker(bind=engine)
    session = Session()
    return engine, session
