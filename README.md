# ncbi_survey
Retrieves genome assembly metadata from the NCBI assembly repository for a list of accessioned assemblies.

***

    usage: VEuPath_taxon_get.py [-h] [--out OUT] [--test] [--assemblies ASSEMBLIES] [--sampleN SAMPLEN]

    Generate report of VEuPath and NCBI/Refseq assembly overlaps.

    options:
        -h, --help            show this help message and exit
        --out OUT             path of output report html file (default=refstatus_report.html)
        --test                run in test mode using internal data
        --assemblies ASSEMBLIES
                        path of assembly input file
        --sampleN SAMPLEN     number of retrieved assemblies to process

***
This program compares [NCBI format assembly accessions](https://www.ncbi.nlm.nih.gov/assembly/help/#accessions) (GCA_*  AND GCF_*) to a supplied list of accessions and
locates all other assemblies available for the linked NCBI taxon id. The program then compares the retrieved
assembly accessions to determine if the supplied accession is either behind the current assembly version, or
is does not match the Refseq nominated reference assembly.

By default the list of NCBI assembly accessions is derived from the VEuPathDB web service 

<https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/standard>

A user supplied list of assembly accessions can be provided with the --assemblies option. The file supplied
to the --assemblies option must have 5 fields delimited with the "|" symbol and conform to the format

`<NCBI assembly accession>|<name>|<is reference yes|no>|<assembly version>|<NCBI taxon id>`

e.g.

GCA_001949145.1|Acanthaster planci (Crown-of-thorns starfish)|yes|OKI_Apl_1.0|133434

Example input files are supplied in the test_data directory. An small inbuilt set of test organisms can
be used with the --test option.

A subset of the supplied assemblies can be processed using the --sampleN <int> option.

For each assembly record the supplied NCBI taxon ID is used to perform an NCBI assembly lookup to obtain all 
other assemblies linked to that taxon id, or it's parent (in the case of strain level assemblies).Aallocation 
of strain level taxon ids is now deprecated at NCBI but the lookup is required when dealing with lists of 
older assemblies).

The output report is generated as an HTML table listing all linked assemblies, along with colour coding to
indicate at what level the assembly accession in the supplied list can be matched to the available records
(assemblies can be updated or deprecated over time). Matching to source (GCA_ or GCF_), assembly id, and
assembly version are all supported.

***

Example usages:

1. Retrieve all VEuPathDB assemblies, default HTML file output

`python3 VEuPath_taxon_get.py`

2. Retrieve assemblies using a supplied list to the output file foo.html

`python3 VEuPath_taxon_get.py --out foo.html --assemblies test_data/tst_ensembl_metazoa_assemblies`

3. Retrieve only the firts 10 assemblies using a supplied list

`python3 VEuPath_taxon_get.py --sampleN 10 --assemblies test_data/tst_ensembl_metazoa_assemblies`   
