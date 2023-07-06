# ncbi_survey
Retrieves genome assembly metadata from the NCBI assembly repository for a list of accessioned assemblies.

usage: VEuPath_taxon_get.py [-h] [--out OUT] [--test] [--sampleN SAMPLEN]

Generate report of VEuPath and NCBI/Refseq assembly overlaps.

options:
  -h, --help         show this help message and exit
  --out OUT          path of output report html file (default={default_output})
  --test             run in test mode using internal data
  --sampleN SAMPLEN  number of retrieved assemblies to process

The list of NCBI assembly accessions is derived from the VEuPathDB web service 

https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/standard

For each assembly record the supplied NCBI taxon ID is used to perform an NCBI assembly lookup to obtain all 
other assemblies linked to that taxon id or it's parent in the case of strain level assemblies (allocation of
strain level taxon ids is now deprecated at NCBI but the lookup is required when dealing with lists of older
assemblies).

The output report is generated as an HTML table listing all linked assemblies, along with colour coding to
indicate at what level the assembly accession in the supplied list can be matched to the available records
(assemblies can be updated or deprecated over time).