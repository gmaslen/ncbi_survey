import typing
import requests
import json
import sys
import re
import argparse
from jinja2 import FileSystemLoader, Environment
from ncbi.datasets.metadata.genome import print_assembly_metadata_by_fields
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_asm_accessions
from ncbi.datasets.metadata.genome import assembly_values_by_fields
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_taxon

class ncbi_accession:
    """generates an object based on an INSDC format assembly accession string. 
    Can be used for equivalence comparisons. Provides accessor methods for substring components
    n.accession - supplied accession (string)
    n.source - prefix code for INSDC source (GCA/GCF) (sting)
    n.id - assembly id (string)
    n.version - assembly version (integer)
    """    
    def __init__(self, fullAcc: str) -> None:
        acc = fullAcc.split('.')
        parts = acc[0].split('_')
        self.source = parts[0]
        self.id = parts[1]
        self.version = acc[1]
        
    def __eq__(self,other) -> bool:
        return (self.source, self.id, self.version) == (other.source, other.id, other.version)
        
    def __ne__(self,other) -> bool:
        return (self.source, self.id, self. version) != (other.source, other.id, other.version)
        
    def accession(self) -> str:
        return str(self.acc)
    def source(self) -> str:
        return str(self.source)
    def id(self) -> str:
        return str(self.id)
    def version(self) -> int:
        return (self.version)

class ncbi_assembly:
    """object storing NCBI assembly properties - available object methods
    self.accession() = retrieve NCBI accessiom
    self.category()  = retrieve NCBI assembly category
    self.source()    = retrieve NCBI assembly source
    self.submitted() = retrieve NCBI assembly date of submission
    self.taxid()     = retrieve NCBI taxonomy id
    self.strain      = retrieve NCBI strain name or None
    self.busco       = refrence to busco class object
    """    
    def __init__(self, accession , category, level, genes, proteins, source, submitted, tax_id, strain,busco ) -> None:
        self.accession = accession
        self.category = category
        self.level = level
        self.genes = genes
        self.proteins = proteins
        self.source = source
        self.submitted = submitted
        self.taxid = tax_id
        self.strain = strain
        self.busco = busco
        
    def accession(self) -> str :
        return self.accession 
    def category(self) -> str :
        return self.category
    def level(self) -> str :
        return self.level
    def genes(self) -> str :
        return self.genes
    def proteins(self) -> str :
        return self.proteins
    def source(self) -> str:
        return self.source
    def submitted(self)  -> str :       
        return self.submitted
    def taxid(self) -> str :
        return self.taxid
    def strain(self) -> str :
        return self.strain
    def busco(self) -> str:
        return self.busco
    
class busco:
    """handler class for NCBI BUSCO statistics
    NCBI BSCO scores are presented as floats not absolute counts (as is normal fro reporting)
    """    
    def __init__(self, stats ) -> None: 
        self.completed = round(stats.get('complete'), 4)
        self.single_copy = round(stats.get('single_copy'), 4)
        self.duplicated = round(stats.get('duplicated'), 4)
        self.fragmented = round(stats.get('fragmented'), 4)
        self.missing = round(stats.get('missing'), 4)
        self.db = stats.get('busco_lineage')
        self.version = stats.get('busco_ver')
        self.total = stats.get('total_count')
                 
    def info(self) -> dict:
        return self
    def db(self) -> str:
        return self.db
    def version(self) -> str:
        return self.busco
    def total(self) -> str:
        return self.total
    def completed(self) -> str:
        return self.completed
    def single_copy(self) -> str:
        return self.single_copy
    def duplicated(self) -> str:
        return self.duplicateded
    def fragmented(self) -> str:
        return self.fragmented
    def missing(self) -> str:
        return self.missing
    def __str__(self):
        return f"{self.db} v{self.version} C:{self.completed}[S:{self.single_copy},D:{self.duplicated}],F:{self.fragmented},M:{self.missing},n:{self.total}"

class VPassembly:
    """generate object to store assembly details retrieved from VEuPath web service.
    Available methods:
    accession = return VEuPath assembly accession
    name      = return VEuPath assembly name
    reference = refturn VEuPath assembly reference genome status
    version   = return VEuPath assembly version
    taxid     = return VEuPath assembly NCBI taxon id
    """    
    def __init__(self, accession , name, reference, version, tax_id ) -> None:
        self.accession = accession
        self.name = name
        self.reference = reference
        self.version = version
        self.taxid = tax_id
    
    def __str__(self):
        return f"name={self.name}|accession={self.accession}|reference={self.reference}|version={self.version}|taxid={self.taxid}"    
    
    def accession(self) -> str :
        return self.accession 
    def name(self) -> str :
        return self.name
    def reference(self) -> str:
        return self.ref
    def version(self)  -> str :       
        return self.version
    def taxid(self) -> str :
        return self.taxid
     
def getVEuPathAssemblies() -> list:
    """Retrieve VEuPathDB assembly information for all organisms using Webservices API

    Returns:
        dict: dictionary of VEuPathDB assembly objects using INSDC assembly accession as key (i.e. GCA_<nnnnn>)
    """
        
    base_url ="https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/standard"
    hm = base_url + '?reportConfig={"attributes":["primary_key","genome_source","genome_version","annotation_version","is_reference_strain","ncbi_tax_id"],"tables":[],"attributeFormat":"text"}'

    x = requests.get(hm)
    assemblies = []

    if x.ok :
        data= x.json()
        taxids = 0
        rows = 0
    
        for l in data.get('records'):
            rows += 1
            v = VPassembly(
                str(l.get('attributes').get('genome_version')),
                str(l.get('attributes').get('primary_key')),
                str(l.get('attributes').get('is_reference_strain') ),
                str(l.get('attributes').get('annotation_version')),
                str(l.get('attributes').get('ncbi_tax_id'))
            )
            assemblies.append(v)
        return assemblies
    else:   
        print( str(x.status_code) )
        print( x.apparent_encoding )
        req = x.request
        print(req.body)

def getNCBIrecord( VPassembly: type [VPassembly] ) -> list :
    """retrieve NCBI assembly record based on assembly accession, then use this to perform a lookup 
    of the linked NCBItaxon id and retrieve all possible assemblies for that taxon. This required as 
    the supplied assembly accession may be linked to a strain specific taxon id

    Args:
        VPassembly (list): List of VPassembly objects containing INSDC GC[AF] accessions to retrieve
    """
    ret = []
    GCAid = [ VPassembly.accession]
    for assembly in get_assembly_metadata_by_asm_accessions( GCAid ):
        NCBIgca = assembly.get('assembly').get('assembly_accession')
        
        vp_ref = VPassembly.reference
        ncbi_ref = assembly.get('assembly').get('assembly_category')
        ncbi_tax_id = 'unknown' if not assembly.get('assembly').get('org').get('tax_id') else assembly.get('assembly').get('org').get('tax_id')
        
        org_assemblies = getAssembliesByTaxon(ncbi_tax_id)
        arrayAnalysis = analyseAssemblies ( org_assemblies, NCBIgca)
        
        fields =  { 
                   "num_assemblies": str( len( org_assemblies) ) , 
                   "ncbi_accession": NCBIgca ,
                   "vp_name": VPassembly.name , 
                   "vp_version": VPassembly.version, 
                   "vp_taxid": VPassembly.taxid,
                   "ncbi_taxid=": ncbi_tax_id
                   }

        for x in arrayAnalysis:
            fcopy = fields.copy()
            fcopy.update(x)
            ret.append(fcopy)

    return ret                
            
def analyseAssemblies( org_assemblies : list, NCBIgca : str) -> list:
    """checks for matches to supplied assembly accession in supplied array of 
    NCBI assembly ids for a taxon. Returns a list of list objects, each of which
    represents an individual NCBI assembly and an analysis of the type of match
    for that assembly. Possible values are.
    
    NCBI non-reference assembly      -  NCBI assembly is non-reference and does 
                                        not match supplied assembly accession
    matched reference assembly       -  NCBI assembly is a reference genome and
                                        matches the supplied assembly accession
    matched non-reference assembly   -  NCBI assembly is non-reference and matches
                                        the supplied assembly accession
    
    Args:
        org_assemblies (list): list of NCBI assembly objects
        NCBIgca (str): NCBI assembly string to match

    Returns:
        list:   list of dicts. Each list object summarises the NCBI assembly
                and whether it matches the supplied assembly accession string
    """    
    retdata =[]
    
    for n in org_assemblies:
            
        cat = "" if not n.category else n.category 
        report_string = "NCBI non-reference assembly"
        vpAcc = ncbi_accession(NCBIgca)
        ncbiAcc = ncbi_accession(n.accession)
        report = checkAccessions( vpAcc , ncbiAcc)
            
        if cat == "representative genome":
            if len(report) > 0 :
                report_string = ":".join(report)
            else:
                report_string = 'matched reference assembly'
        elif cat == "" :
            if len(report) == 0:
                report_string = 'matched non-reference assembly'
    
        r= {
            "ncbi_accession": str(n.accession), 
            "ncbi_category": str(n.category),
            "ncbi_level": str(n.level),
            "ncbi_genes": str(n.genes),
            "ncbi_proteins": str(n.proteins), 
            "ncbi_source": str(n.source), 
            "ncbi_submitted": str(n.submitted), 
            "ncbi_taxid": str(n.taxid),
            "ncbi_strain": str(n.strain),
            "ncbi_match": report_string,
            "ncbi_busco": n.busco
        }

        retdata.append(r)
    return retdata

def getAssembliesByTaxon (ncbi_tax_id : str ) -> list:
    """Retrieve all NCBI assemblies for a given NCBI taxon id

    Args:
        ncbi_tax_id (str): NCBI taxon id to use for assembly list retieval

    Returns:
        list: list of ncbi_assembly objects for input taxon id
    """    
    ncbi_obj = []
    for assembly in get_assembly_metadata_by_taxon( ncbi_tax_id ):
        ass_accession = assembly.get('assembly').get('assembly_accession')
        ass_category = assembly.get('assembly').get('assembly_category')
        ass_level = assembly.get('assembly').get('assembly_level')
        submitter = assembly.get('assembly').get('submitter')
        submission_date = assembly.get('assembly').get('submission_date')
        orgdat = getOrgData( assembly )
        ncbi_tax_id = orgdat.get('tax_id')
        strain = orgdat.get('strain')
        annotation = getAnnoMetadata( assembly )
        total_genes = annotation.get('total_genes') if annotation.get('total_genes') else 0
        total_proteins = annotation.get('total_proteins') if annotation.get('total_proteins') else 0
        busco_obj = getBuscoStats( assembly )
        if busco_obj :
            print( str(ass_accession) + " has BUSCO scores " + str(busco_obj) )
        ncbi_obj.append( ncbi_assembly( ass_accession, ass_category, ass_level, total_genes, total_proteins, submitter, submission_date, ncbi_tax_id, strain, busco_obj) )
        
    return ncbi_obj

def getAnnoMetadata( assembly:dict ) -> dict :
    """parse required assembly annotation data from NCBI assembly dict

    Args:
        assembly (dict): NCBI assembly dictionary

    Returns:
        dict: fields extracted from NCBI annotation_metadata node in the NCBI dictionary
    """  
    ret ={}
    if assembly.get('assembly').get('annotation_metadata'):
        annstats = assembly.get('assembly').get('annotation_metadata').get('stats').get('gene_counts')
        ret.update({ 'total_genes': annstats.get('total') if annstats.get('total') else 0 })
        ret.update({'total_proteins': annstats.get('protein_coding') if annstats.get('protein_coding') else 0 })
    return ret
        
def getOrgData( assembly:dict) -> dict :
    """parse required organism data from NCBI assembly dict

    Args:
        assembly (dict): NCBI assembly dictionary

    Returns:
        dict: fields extracted from NCBI org node in the NCBI dictionary
    """    
    ret ={}
    if assembly.get('assembly').get('org'):
        orgdat = assembly.get('assembly').get('org')
        ret.update( {'strain': orgdat.get('strain') if orgdat.get('strain') else 'unknown' } )
        ret.update( {'tax_id': orgdat.get('tax_id')} )
        ret.update( {'parent_tax_id': orgdat.get('parent_tax_id')} )
    return ret
    
def getBuscoStats( assembly: dict ) -> dict :
    """extract Busco score values from NCBI assembly object

    Args:
        assembly (dict): NCBI assembly object dictionary

    Returns:
        dict: dictionary of Busco values
    """    
    
    if assembly.get('assembly').get('annotation_metadata'):
        if assembly.get('assembly').get('annotation_metadata').get('busco'):
            bstats = assembly.get('assembly').get('annotation_metadata').get('busco')
            b = busco(bstats)
            return b
        else:
            return None
    else:
        return None

def checkAccessions( a : type [ncbi_accession] , b : type [ncbi_accession] ) -> list:
    """compare two INSDC format assembly accession strings and return list of errors detected

    Args:
        a (type[ ncbi_accession]): ncbi_accession object
        b (type[ ncbi_accession ]): ncbi_accession object

    Returns:
        list: list of mismatch types detected between the two accessions. Empty set equals identity.
    """    
    r=[]
    if ( a != b):
        if( a.source != b.source):
            r.append( 'source mismatch')
        if( a.id != b.id ):
            r.append('id mismatch')
        if( a.version != b.version ):
            r.append('version mismatch')
    return r
            
def read_file( fh_in:str)-> list :
    """Read | delimited input file of assemblies of the format
    
    <NCBI assembly accession> , <name>, <reference yes|no>, <assembly version (int)>, <NCBI taxon id>
    
    e.g.
    GCF_943734845.2|Anopheles funestus|yes|AfunF3.2|62324
    
    Args:
        fh_in (str): input filepath

    Returns:
        list: array of assembly objects
    """
    ret = []
    with open(fh_in, "r") as fh:
        for line in fh:
            l  =line.strip()
            ta = l.rsplit('|')
            
            o = VPassembly( ta[0] , ta[1], ta[2], ta[3], ta[4] )
            ret.append(o)

    return ret
    
    
def main():
    
    parser = argparse.ArgumentParser(description='Generate report of VEuPath and NCBI/Refseq assembly overlaps.')
    default_output = 'refstatus_report.html'
    parser.add_argument('--out', help='path of output report html file (default={default_output})', type=str, required=False )
    parser.add_argument('--test', help=f'run in test mode using internal data', action="store_true" )
    parser.add_argument('--assemblies', help=f'path of assembly input file', type=str, required=False )
    parser.add_argument('--sampleN', help=f'number of retrieved assemblies to process', type=int, required=False )
    args = parser.parse_args()

    if args.out:
        default_output = args.out
    
    if args.test:
        #array of test data for debugging   
        tst_busco = VPassembly( 'GCF_943734845.2', 'Anopheles funestus', 'yes', 'AfunF3.2', '62324' )
        tst_strain = VPassembly( 'GCA_003013265.1', 'Trypanosoma congolense' , 'no', 'unknown', '1068625')
        assemblies = [ tst_busco, tst_strain]
        print("using internal test data with " + str(len(assemblies)) + " items")
    elif args.assemblies :
        assemblies = read_file( args.assemblies )
    else:
        assemblies = getVEuPathAssemblies()
        print( str(len(assemblies)) + " records retrieved from VEuPath")
    
    skipped = [] 
    allrecords = []
    failed = []
    

    
#use assemblies[:n] to restrict processing to n entries while testing
    
    if args.sampleN :
        subsample = args.sampleN
        print( "processing " + str(subsample) + " of " + str(len(assemblies)) + " records" )
    else:
        subsample = len(assemblies)
        
    for a in assemblies[:subsample]:

        if re.search( "GC[AF]_" , str(a.accession) ) :
            r = getNCBIrecord( a )
            if len(r) == 0:
                failed.append(a)
            else:
                allrecords.append(r)
        else:
            print("Invalid genome accession format for entry - skipping - ")
            print( a )
            skipped.append(a)
    
    print("Total skipped records = " + str(len(skipped)) )
    print("Total failed retrieval records = " + str(len(failed)) )
            
    htmlAssemblyReport( default_output, allrecords, failed, skipped)
    
def htmlAssemblyReport( fh_out: str, ASMlist: list, failed: list, skipped: list ):
    """Print HTML table summary of retrieved assembly records

    Args:
        ASMlist (list): list of assembly data hashes retrieved
    """    
    fh = open( fh_out, "w")
    environment = Environment(loader=FileSystemLoader("templates/"))
    template = environment.get_template("base.html")
    context = { "record": ASMlist, "failed_records": failed, "skipped_records": skipped }
    fh.write( template.render(context) )
    fh.close()
           

if __name__ == '__main__':
    main()
    sys.exit(0)