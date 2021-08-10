/*
========================================================================================
    Config file for defining SSDS intervals defined in previous work
========================================================================================

*/
params {
    ssds_intervals {
        mm10 {
            hotspots = [
                        "13R":'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954833&format=file&file=GSM1954833%5F13R%5Fhotspots%2Etab%2Egz',
                        "B6":'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954835&format=file&file=GSM1954835%5FB6%5Fhotspots%2Etab%2Egz',
                        "C3H":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954841&format=file&file=GSM1954841%5FC3H%5Fhotspots%2Etab%2Egz",
                        "CAST":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954846&format=file&file=GSM1954846%5FCAST%5Fhotspots%2Etab%2Egz",
                        "MOL":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954855&format=file&file=GSM1954855%5FMOL%5Fhotspots%2Etab%2Egz",
                        "PWD":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954858&format=file&file=GSM1954858%5FPWD%5Fhotspots%2Etab%2Egz",
                        "B6 PRDM9ko":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1954864&format=file&file=GSM1954864%5FB6%5FPrdm9ko%5Fhotspots%2Etab%2Egz"
                        ]
                        
            otherloci = [
                        "origins":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE148327&format=file&file=GSE148327%5Fhiconf%5Forigins%2Emm10%2Ebedgraph%2Egz"
                        ]
        }
        hg38 {
            hotspots = [
                        "A/A1":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072244&format=file&file=GSM5072244%5FAA1%2Epeaks%2Ebedgraph%2Egz",
                        "A/A2":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072246&format=file&file=GSM5072246%5FAA2%2Epeaks%2Ebedgraph%2Egz",
                        "A/B":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072254&format=file&file=GSM5072254%5FAB1%2Epeaks%2Ebedgraph%2Egz",
                        "A/C":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072256&format=file&file=GSM5072256%5FAC%2Epeaks%2Ebedgraph%2Egz",
                        "A/N":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072258&format=file&file=GSM5072258%5FAN%2Epeaks%2Ebedgraph%2Egz",
                        "C/L4":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5072260&format=file&file=GSM5072260%5FCL4%2Epeaks%2Ebedgraph%2Egz",
                        "ALL":"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166483&format=file&file=GSE166483%5FhotspotsData%2Etab%2Egz"
                        ]
        }
    }
}