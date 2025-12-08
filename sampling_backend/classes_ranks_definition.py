from enum import Enum

'''---taxonomic rank model organisation---'''
class RankType(Enum):
    DOMAIN = "DOMAIN"
    KINGDOM = "KINGDOM"
    SUBKINGDOM = "SUBKINGDOM"
    SUPERPHYLUM = "SUPERPHYLUM"
    PHYLUM = "PHYLUM"
    SUBPHYLUM = "SUBPHYLUM"
    INFRAPHYLUM = "INFRAPHYLUM"
    SUPERCLASS = "SUPERCLASS"
    CLASS = "CLASS"
    SUBCLASS = "SUBCLASS"
    INFRACLASS = "INFRACLASS"
    COHORT = "COHORT"
    SUBCOHORT = "SUBCOHORT"
    SUPERORDER = "SUPERORDER"
    ORDER = "ORDER"
    SUBORDER = "SUBORDER"
    INFRAORDER = "INFRAORDER"
    PARVORDER = "PARVORDER"
    SUPERFAMILY = "SUPERFAMILY"
    FAMILY = "FAMILY"
    SUBFAMILY = "SUBFAMILY"
    TRIBE = "TRIBE"
    SUBTRIBE = "SUBTRIBE"
    GENUS = "GENUS"
    SUBGENUS = "SUBGENUS"
    SECTION = "SECTION"
    SUBSECTION = "SUBSECTION"
    SERIES = "SERIES"
    SUBSERIES = "SUBSERIES"
    SPECIES_GROUP = "SPECIES_GROUP"
    SPECIES_SUBGROUP = "SPECIES_SUBGROUP"
    SPECIES = "SPECIES"
    FORMA_SPECIALIS = "FORMA_SPECIALIS"
    SUBSPECIES = "SUBSPECIES"
    VARIETAS = "VARIETAS"
    SUBVARIETY = "SUBVARIETY"
    FORMA = "FORMA"
    SEROGROUP = "SEROGROUP"
    SEROTYPE = "SEROTYPE"
    STRAIN = "STRAIN"
    ISOLATE = "ISOLATE"
    # Extra ranks from datasets (used as "unranked" in Taxallnomy-style logic)
    CLADE = "CLADE"
    NO_RANK = "NO_RANK"
    SUPERKINGDOM = "SUPERKINGDOM"
    REALM = "REALM"
    BIOTYPE = "BIOTYPE"
    GENOTYPE = "GENOTYPE"
    MORPH = "MORPH"
    PATHOGROUP = "PATHOGROUP"
    ACELLULAR_ROOT = "ACELLULAR_ROOT"
    CELLULAR_ROOT = "CELLULAR_ROOT"

RANK_ORDER = [
    RankType.DOMAIN,
    RankType.KINGDOM,
    RankType.SUBKINGDOM,
    RankType.SUPERPHYLUM,
    RankType.PHYLUM,
    RankType.SUBPHYLUM,
    RankType.INFRAPHYLUM,
    RankType.SUPERCLASS,
    RankType.CLASS,
    RankType.SUBCLASS,
    RankType.INFRACLASS,
    RankType.COHORT,
    RankType.SUBCOHORT,
    RankType.SUPERORDER,
    RankType.ORDER,
    RankType.SUBORDER,
    RankType.INFRAORDER,
    RankType.PARVORDER,
    RankType.SUPERFAMILY,
    RankType.FAMILY,
    RankType.SUBFAMILY,
    RankType.TRIBE,
    RankType.SUBTRIBE,
    RankType.GENUS,
    RankType.SUBGENUS,
    RankType.SECTION,
    RankType.SUBSECTION,
    RankType.SERIES,
    RankType.SUBSERIES,
    RankType.SPECIES_GROUP,
    RankType.SPECIES_SUBGROUP,
    RankType.SPECIES,
    RankType.FORMA_SPECIALIS,
    RankType.SUBSPECIES,
    RankType.VARIETAS,
    RankType.SUBVARIETY,
    RankType.FORMA,
    RankType.SEROGROUP,
    RankType.SEROTYPE,
    RankType.STRAIN,
    RankType.ISOLATE
]

RANK_ABBREVIATIONS = {
    RankType.DOMAIN: "Dom", RankType.KINGDOM: "Kin", RankType.SUBKINGDOM: "sbKin",
    RankType.SUPERPHYLUM: "spPhy", RankType.PHYLUM: "Phy", RankType.SUBPHYLUM: "sbPhy", RankType.INFRAPHYLUM: "inPhy",
    RankType.SUPERCLASS: "spCla", RankType.CLASS: "Cla", RankType.SUBCLASS: "sbCla", RankType.INFRACLASS: "inCla",
    RankType.COHORT: "Coh", RankType.SUBCOHORT: "sbCoh",
    RankType.SUPERORDER: "spOrd", RankType.ORDER: "Ord", RankType.SUBORDER: "sbOrd",
    RankType.INFRAORDER: "inOrd", RankType.PARVORDER: "prOrd",
    RankType.SUPERFAMILY: "spFam", RankType.FAMILY: "Fam", RankType.SUBFAMILY: "sbFam",
    RankType.TRIBE: "Tri", RankType.SUBTRIBE: "sbTri",
    RankType.GENUS: "Gen", RankType.SUBGENUS: "sbGen",
    RankType.SECTION: "Sec", RankType.SUBSECTION: "sbSec",
    RankType.SERIES: "Ser", RankType.SUBSERIES: "sbSer",
    RankType.SPECIES_GROUP: "Sgr", RankType.SPECIES_SUBGROUP: "sbSgr",
    RankType.SPECIES: "Spe", RankType.FORMA_SPECIALIS: "Fsp",
    RankType.SUBSPECIES: "sbSpe", RankType.VARIETAS: "Var",
    RankType.SUBVARIETY: "sbVar", RankType.FORMA: "For",
    RankType.SEROGROUP: "Srg", RankType.SEROTYPE: "Srt",
    RankType.STRAIN: "Str", RankType.ISOLATE: "Iso"
}

UNRANKED_TYPES = {
    RankType.CLADE,
    RankType.NO_RANK,
    RankType.BIOTYPE,
    RankType.GENOTYPE,
    RankType.MORPH,
    RankType.PATHOGROUP,
    RankType.ACELLULAR_ROOT,
    RankType.CELLULAR_ROOT
}

'''---core classes---'''
class TaxonNode:
    def __init__(self, tax_id, parent_id, original_rank, rank, name, phantom=False, interpolated=False):
        self.tax_id = tax_id
        self.parent_id = parent_id
        self.original_rank = original_rank
        self.rank = rank
        self.rank_level = RANK_ORDER.index(rank) if rank in RANK_ORDER else None
        self.name = name
        self.original_name = name  # for user lookups / pre-interpolation name

        self.sibling = None
        self.children = []
        self.genomes = []

        self.span = (None, None)  # leaf index span (from linearize_tree)
        self.genome_span = (None, None)  # genome index span (from annotate_tree)
        self.nearest_descendant_rank_level = None
        self.nearest_descendant_rank = None

        self.phantom = phantom
        self.interpolated = interpolated

        # anomaly flags
        self.is_self_parent = False
        self.is_orphan = False
        self.is_unreachable = False

        # summaries attached to root; None for most nodes
        self.anomalies_total = None
        self.genome_source_stats = None

    def add_child(self, child_node):
        self.children.append(child_node)

    def add_genome(self, genome_record):
        self.genomes.append(genome_record)

    def genome_count(self):
        start, end = getattr(self, "genome_span", (None, None))
        if start is not None and end is not None:
            return end - start
        return 0

class GenomeRecord:
    def __init__(self, accession, tax_id, name, refseq_category, assembly_level, is_reference, origin_source=None,
        source_database=None,
        ncbi_source_database=None,
        ensembl_division=None,
        assembly_name=None,
        ensembl_core_db=None,):
        self.accession = accession
        self.tax_id = tax_id
        self.name = name
        self.refseq_category = refseq_category
        self.assembly_level = assembly_level
        self.is_reference = is_reference
        self.rank = None
        self.rank_level = None

        # source info
        self.origin_source = origin_source  # "NCBI" or "Ensembl" (or None)
        self.source_database = source_database  # e.g. "SOURCE_DATABASE_REFSEQ"

        # Ensembl-specific data
        self.ensembl_division = ensembl_division  # e.g. "EnsemblVertebrates"
        self.assembly_name = assembly_name  # Ensembl assembly name string
        self.ensembl_core_db = ensembl_core_db  # core DB name



