import csv
import requests
import logging
import sys
import re
import argparse
from jinja2 import FileSystemLoader, Environment
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_asm_accessions
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_taxon

## unusued imports:
# from ncbi.datasets.metadata.genome import print_assembly_metadata_by_fields
# from ncbi.datasets.metadata.genome import assembly_values_by_fields

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
        self.complete = stats.get('complete', 0.0)
        self.single_copy = stats.get('single_copy', 0.0)
        self.duplicated = stats.get('duplicated', 0.0)
        self.fragmented = stats.get('fragmented', 0.0)
        self.missing = stats.get('missing', 0.0)
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
    def complete(self) -> str:
        return self.complete
    def single_copy(self) -> str:
        return self.single_copy
    def duplicated(self) -> str:
        return self.duplicated
    def fragmented(self) -> str:
        return self.fragmented
    def missing(self) -> str:
        return self.missing
    def __str__(self):
        return f"{self.db} v{self.version} C:{self.complete:.3f}[S:{self.single_copy:.3f},D:{self.duplicated:.3f}],F:{self.fragmented:.3f},M:{self.missing:.3f},n:{self.total}"

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

    asm_objs = requests.get(hm)
    assemblies = []

    if asm_objs.ok:
        data= asm_objs.json()
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
        print( str(asm_objs.status_code) )
        print( asm_objs.apparent_encoding )
        req = asm_objs.request
        print(req.body)

def getNCBIrecord( Original_Accn: str, VPassembly: VPassembly ) -> list :
    """Retrieve NCBI assembly record based on assembly accession, then use this to perform a lookup 
    of the linked NCBItaxon id and retrieve all possible assemblies for that taxon. This required as 
    the supplied assembly accession may be linked to a strain specific taxon id

    Args:
        VPassembly: List of VPassembly objects containing INSDC GC[AF] accessions to retrieve
    """
    ret = []
    GCAid = [ VPassembly.accession ]

    for assembly in get_assembly_metadata_by_asm_accessions( GCAid ):
        # print(assembly)
        NCBIgca = assembly.get('assembly').get('assembly_accession')

        ncbi_tax_id = 'unknown' if not assembly.get('assembly').get('org').get('tax_id') else assembly.get('assembly').get('org').get('tax_id')

        org_assemblies, mismatch_version_meta = getAssembliesByTaxon(Original_Accn, ncbi_tax_id)
        arrayAnalysis = analyseAssemblies ( org_assemblies, NCBIgca)

        fields =  {
                   "num_assemblies": str( len( org_assemblies) ),
                   "ncbi_accession": NCBIgca,
                   "vp_name": VPassembly.name,
                   "vp_version": VPassembly.version,
                   "vp_taxid": VPassembly.taxid,
                   "ncbi_taxid": ncbi_tax_id
                   }

        for x in arrayAnalysis:
            fcopy = fields.copy()
            fcopy.update(x)
            ret.append(fcopy)

    return ret, mismatch_version_meta             
            
def analyseAssemblies( org_assemblies : list, NCBIgca : str) -> list:
    """Checks for matches to supplied assembly accession in supplied array of 
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
        org_assemblies: list of NCBI assembly objects
        NCBIgca: NCBI assembly string to match

    Returns:
        list of dicts. Each list object summarises the NCBI assembly
        and whether it matches the supplied assembly accession string
    """    

    retdata =[]
    for n in org_assemblies:
        
        ## Checks if assembly has representative/reference genome status
        cat = "" if not n.category else n.category

        report_string = "NCBI non-reference assembly"
        vpAcc = ncbi_accession(NCBIgca)
        ncbiAcc = ncbi_accession(n.accession)
        report = checkAccessions(vpAcc, ncbiAcc)
            
        if cat == "representative genome":
            if len(report) > 0 :
                report_string = "Acc mismatch [" + ", ".join(report) + "]"
            else:
                report_string = 'User Acc matched reference assembly'
        elif cat == "":
            cat = "non-reference"
            if len(report) == 0:
                report_string = 'User Acc matched non-reference assembly'
    
        r= {
            "ncbi_accession": str(n.accession), 
            "ncbi_category": cat.title(),
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

def getAssembliesByTaxon (Original_Accn: str, ncbi_tax_id: str ) -> list:
    """Retrieve all NCBI assemblies for a given NCBI taxon id

    Args:
        ncbi_tax_id: NCBI taxon id to use for assembly list retieval

    Returns:
        List of ncbi_assembly objects for input taxon id
    """    
    ncbi_obj = []
    meta_version_mismatch = None

    for assembly in get_assembly_metadata_by_taxon( ncbi_tax_id ):
        asm_accession = assembly.get('assembly').get('assembly_accession')
        asm_category = assembly.get('assembly').get('assembly_category')
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
        if busco_obj:
            print( str(asm_accession) + " has BUSCO scores " + str(busco_obj) )
        ncbi_obj.append( ncbi_assembly( asm_accession, asm_category, ass_level, total_genes, total_proteins, submitter, submission_date, ncbi_tax_id, strain, busco_obj) )

        # Testing whether or not api call has pulled back the query accession meta or a later asm version
        if asm_accession.startswith(Original_Accn.split(".")[0]) and not asm_accession.endswith(Original_Accn.split(".")[1]):
            logging.warning(f"Query accession {Original_Accn} not returned via API ! [Version obtained: {asm_accession}]")
            meta_version_mismatch = [Original_Accn, asm_accession]
    
    return ncbi_obj, meta_version_mismatch

def getAnnoMetadata( assembly:dict ) -> dict :
    """parse required assembly annotation data from NCBI assembly dict

    Args:
        assembly: NCBI assembly dictionary

    Returns:
        Fields extracted from NCBI annotation_metadata node in the NCBI dictionary
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

def checkAccessions( user_supplied_acc : ncbi_accession , ncbi_acc : ncbi_accession ) -> list:
    """compare two INSDC format assembly accession strings and return list of errors detected

    Args:
        user_supplied_acc: User supplied assembly accession (ncbi_accession object)
        ncbi_acc: Taxon ID derived alternate assembly accession (ncbi_accession object)

    Returns:
        List of mismatch types detected between the two assembly accessions. Empty set equals accession identity.
    """    
    accession_mismatch=[]

    if ( user_supplied_acc != ncbi_acc):
        if( user_supplied_acc.source != ncbi_acc.source):
            accession_mismatch.append('Source') # GenBank Vs RefSeq mismatch
        if( user_supplied_acc.id != ncbi_acc.id ):
            accession_mismatch.append('ID') # Assembly GCA/GFC ID GC{CF}_ id -> 'XXXXXXXXX'
        if( user_supplied_acc.version != ncbi_acc.version ):
            accession_mismatch.append('Version') # Assembly version mismatch
    return accession_mismatch
            
def read_assembly_queries(fh_in:str)-> list :
    """Read | delimited input file of assemblies of the format
    
    <NCBI assembly accession> , <name>, <reference yes|no>, <assembly version (int)>, <NCBI taxon id>
    
    e.g.
    GCF_943734845.2|Anopheles funestus|yes|AfunF3.2|62324
    
    Args:
        fh_in: input filepath

    Returns:
        lList of array of assembly objects
    """
    ret = []
    with open(fh_in, "r") as fh:
        for line in fh:
            l  =line.strip()
            ta = l.rsplit('|')
            
            o = VPassembly( ta[0] , ta[1], ta[2], ta[3], ta[4] )
            ret.append(o)

    return ret

def tsvAssemblyReport( output_file: str, asm_records: list, match_type: bool):
    """A function to iterate over all assembly match records returned from a successfully
    user accession lookup. Matches include user supplied assembly accession hits and also
    non-user qeury asm accns identified via NCBI taxon_id lookup.

    Args:
        output_file: Default/user specified output file prefix
        asm_records: Set of all asm record lookups
        match_type: Bool setter for result outfile suffix, True = user accession matched
    """

    header_list = ['Total assemblies','NCBI accession','VEuPath Name','VEuPath version',
                   'VEuPath Taxon ID','NCBI Taxon ID','NCBI Asm Status','Assembly level',
                   'NCBI Total Genes','NCBI Total Proteins','NCBI Source','Submission date',
                   'NCBI Strain','NCBI Match','BUSCO']

    # Define appropriate outfile suffix depending on match_type result set
    if match_type:
        query_acc_output_tsv = f"{output_file}.TargetAccnQueries.tsv"
    else:
        query_acc_output_tsv = f"{output_file}.OtherAssociated-asms-on-taxid.tsv"

    with open(query_acc_output_tsv, 'w+') as tsv_out:
        writer = csv.writer(tsv_out, delimiter='\t', lineterminator='\n')
        writer.writerow(header_list)

        #Iterate over all assemblies (>=1) result sets & write individual accession result tsv
        for asm_matched in asm_records:
            writer.writerow(asm_matched.values())

        tsv_out.close()

def categorize_results (asm_set: list) -> list:
    """Function to separate out sets of results, based on whether or not a specific assembly
    accession was input by the user. Matched set (user query accessions) and
    unmatched (other associated taxon asm accessions)
    
    Args:
        asm_set: Set of assembly result sets passed as [list of {dicts,..}]
    
    Returns:
        matched: Set of asm query accessions matched to reference assembly
        unmatched: Set of asm query accessions not matched to reference
        set_of_references: List of user query accession(s) and associated assembly hit
    """
    matched = []
    unmatched = []
    set_of_references = []

    #Iterate over all assemblies (>=1) result sets & write individual accession result tsv
    for master_asm_index in range(len(asm_set)):
        for asm in asm_set[master_asm_index]:

            match_type = asm.get('ncbi_match')
            ncbi_match = re.search("[Uu]ser Acc", match_type)
            species = asm.get('vp_name')
            asm_accn = asm.get('ncbi_accession')
            ncbi_status = asm.get('ncbi_category')
            reference_or_not = re.search("Representative Genome", ncbi_status)

            if reference_or_not and reference_or_not.group():
                have_reference = True
            else:
                have_reference = False

            if ncbi_match and have_reference == True:
                reference = [species, asm_accn, "Accn matches current reference"]
                set_of_references.append(reference)
                matched.append(asm)

            if ncbi_match and have_reference == False:
                reference = [species, asm_accn, "Non-reference"]
                set_of_references.append(reference)
                matched.append(asm)
            else:
                unmatched.append(asm)

    return matched, unmatched, set_of_references

def htmlAssemblyReport( Output_file: str, ASMlist: list, Failed: list, Skipped: list ):
    """Print HTML table summary of retrieved assembly records based on set of input
    user defined assembly accessions (INSDC accessions i.e. GCA_/GCF_)

    Args:
        Output_file: Name of output file (default: 'refstatus_report.html')
        ASMlist: List of assembly records retrieved from NCBI
        Failed: List of assemblies without NCBI records based on provided accession. 
        Skipped: List of assemblies with improperly defined INSDC accession.
    """    
    fh = open( Output_file, "w")
    environment = Environment(loader=FileSystemLoader("templates/"))
    template = environment.get_template("base.html")
    context = { "record": ASMlist, "failed_records": Failed, "skipped_records": Skipped }
    fh.write( template.render(context) )
    fh.close()

## Main runnable starts here
if __name__ == '__main__':

    skipped = [] 
    allrecords = []
    failed = []
    mismatched_meta = []
    valid_accessions={}
    default_output = 'refstatus_report.html'
    
    parser = argparse.ArgumentParser(description='Generate report of VEuPath and NCBI/Refseq assembly overlaps.')
    parser.add_argument('--out', help='Name of output report html file.', default=default_output, type=str, required=False )
    parser.add_argument('--test', help=f'run in test mode using internal data', action="store_true" )
    parser.add_argument('--assemblies', help=f'path of assembly input file', type=str, required=False )
    parser.add_argument('--sampleN', help=f'number of retrieved assemblies to process', type=int, required=False )
    args = parser.parse_args()

    # Check proper formatting of user output file
    if args.out:
        suffix_match = re.search(".html", args.out)
        if suffix_match:
            default_output = args.out
            file_prefix = default_output.split(".")
        else:
            default_output = args.out + ".html"
            file_prefix = args.out

    if args.test:
        #array of test data for debugging   
        tst_busco = VPassembly( 'GCF_943734845.2', 'Anopheles funestus', 'yes', 'AfunF3.2', '62324' )
        tst_strain = VPassembly( 'GCA_003013265.1', 'Trypanosoma congolense' , 'no', 'unknown', '1068625')
        assemblies = [ tst_busco, tst_strain]
        print("Using internal test data with " + str(len(assemblies)) + " items")
    elif args.assemblies :
        assemblies = read_assembly_queries( args.assemblies )
        print( "Discovered " + str(len(assemblies)) + " user assembly queries in file " + args.assemblies + "\n\n")
    else:
        assemblies = getVEuPathAssemblies()
        print( str(len(assemblies)) + " records retrieved from VEuPath!")

    #use assemblies[:n] to restrict processing to n entries while testing
    if args.sampleN :
        subsample = args.sampleN
        print( "Processing with sub-sample enabled. Operating over " + str(subsample) + " of " + str(len(assemblies)) + " original input assembly queries from file " + str(args.assemblies) )
    else:
        subsample = len(assemblies)

    ## Parse and store information on valid GCA/GCF accessions
    for assembly in assemblies[:subsample]:

        if re.search( "GC[AF]_" , str(assembly.accession) ):
            accession = str(assembly.accession)
            valid_accessions[accession] = assembly
        else:
            print(f"Invalid genome accession format for entry - skipping - {assembly}")
            skipped.append(assembly)

    # Overwrite file of mismathed asm version as as reult of api return
    with open(f"{file_prefix}.api_return_version_mismatched.tsv", 'w') as tsv_out:
        tsv_out.write("#Species\tQuery Accn\tCurrent NCBI Version\n")    
    tsv_out.close
    
    ## Iterate on valid asm accessions
    for valid_assembly_accn, assembly_vp_object in valid_accessions.items():
        species_name = str(assembly_vp_object.name)
        print(f"\nProcessing: [{species_name}] - Accn [{valid_assembly_accn}]")

        asm_record, mismatched_meta = getNCBIrecord(valid_assembly_accn, assembly_vp_object)

        if len(asm_record) == 0:
            failed.append(assembly_vp_object)
        else:
            allrecords.append(asm_record)

        # check for mismatched meta on api call
        if mismatched_meta is not None:
            mismatched_meta.insert(0,f"{species_name}")
            with open(f"{file_prefix}.api_return_version_mismatched.tsv", 'a') as tsv_out:
                writer = csv.writer(tsv_out, delimiter='\t', lineterminator='\n')
                writer.writerow(mismatched_meta)
            tsv_out.close

            
    # Split assembly result sets to isolate user query accessions
    user_matched, user_unmatched, reference_asm_set = categorize_results(allrecords)

    # Summarize output across all user defined accessions with respect to their representative status
    with open(f"{file_prefix}.Asm-Reference-status.tsv", 'w') as tsv_out:
        tsv_out.write("#Species\tQuery Accn\tReference Status\n")
        
        for user_assemblies in reference_asm_set:
            writer = csv.writer(tsv_out, delimiter='\t', lineterminator='\n')
            writer.writerow(user_assemblies)
    tsv_out.close

    # Write output of analysis to HTML output file
    htmlAssemblyReport(default_output, allrecords, failed, skipped)
    
    # For both asm result sets write to TSV    
    tsvAssemblyReport(file_prefix, user_matched, True)
    tsvAssemblyReport(file_prefix, user_unmatched, False)

    print("\nTotal skipped records = " + str(len(skipped)) )
    print("Total failed retrieval records = " + str(len(failed)) )
    
sys.exit(0)